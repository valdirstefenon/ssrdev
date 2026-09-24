import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from Bio import SeqIO
import re
import pandas as pd
from io import StringIO
import os
import numpy as np
from packaging import version
from Bio import __version__ as bio_version

PRIMER3_AVAILABLE = False
try:
    import primer3
    PRIMER3_AVAILABLE = True
except ImportError:
    try:
        from primer3 import primer3 as primer3_module
        primer3 = primer3_module
        PRIMER3_AVAILABLE = True
    except ImportError:
        primer3 = None
        print("WARNING: primer3 module not available. Install with: pip install primer3-py")

from Bio.Blast import NCBIWWW, NCBIXML
import time
import threading
import queue
import logging
import json
import concurrent.futures
import openpyxl
from openpyxl.utils import get_column_letter
from openpyxl.styles import Font
from openpyxl.chart import BarChart, Reference
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
from sklearn.metrics import pairwise_distances
from collections import defaultdict
import gffutils
import gtfparse
import socket

IUPAC_COMPLEMENT = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'N': 'N', 'R': 'Y', 'Y': 'R', 'W': 'W',
    'S': 'S', 'M': 'K', 'K': 'M', 'B': 'V',
    'V': 'B', 'D': 'H', 'H': 'D',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
    'n': 'n', 'r': 'y', 'y': 'r', 'w': 'w',
    's': 's', 'm': 'k', 'k': 'm', 'b': 'v',
    'v': 'b', 'd': 'h', 'h': 'd',
}

class SSRdevApp:
    def __init__(self, root):
        self.root = root
        self.root.title("SSRdev 4.5")
        self.root.geometry("1200x800")

        if not PRIMER3_AVAILABLE:
            messagebox.showwarning(
                "Warning",
                "The 'primer3' module is not installed.\n\n"
                "Primer design will be disabled.\n"
                "To enable primer design, install:\n"
                "python -m pip install primer3-py\n\n"
                "Continuing without primer design..."
            )
            logging.warning("primer3 module not available - primer design disabled")

        if version.parse(bio_version) < version.parse("1.79"):
            messagebox.showerror(
                "Error",
                "This script requires Biopython version 1.79 or higher"
            )
            self.root.destroy()
            return

        # Configurar timeout para conexões de rede
        socket.setdefaulttimeout(30)

        log_path = os.path.join(
            os.getcwd(),
            f"ssrdev_{time.strftime('%Y%m%d_%H%M%S')}.log"
        )
        logging.basicConfig(
            filename=log_path,
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        logging.info("SSRdev 4.5 started")
        if not PRIMER3_AVAILABLE:
            logging.warning("primer3 module not available - primer design disabled")

        self.analysis_in_progress = False

        self._gui_queue: queue.Queue = queue.Queue()
        self._poll_gui_queue()

        self.min_repeats = {
            'mono':  {'motif_length': 1, 'min_repeats': 10, 'enabled': True},
            'di':    {'motif_length': 2, 'min_repeats': 5,  'enabled': True},
            'tri':   {'motif_length': 3, 'min_repeats': 4,  'enabled': True},
            'tetra': {'motif_length': 4, 'min_repeats': 3,  'enabled': True},
            'penta': {'motif_length': 5, 'min_repeats': 3,  'enabled': True},
            'hexa':  {'motif_length': 6, 'min_repeats': 3,  'enabled': True},
        }

        self._ssr_patterns: dict = {}
        self._compile_ssr_patterns()

        self.primer_params = {
            'PRIMER_OPT_SIZE':                    20,
            'PRIMER_MIN_SIZE':                    15,
            'PRIMER_MAX_SIZE':                    30,
            'PRIMER_OPT_TM':                      60.0,
            'PRIMER_MIN_TM':                      50.0,
            'PRIMER_MAX_TM':                      70.0,
            'PRIMER_MIN_GC':                      30.0,
            'PRIMER_MAX_GC':                      70.0,
            'PRIMER_OPT_GC_PERCENT':              50.0,
            'PRIMER_MAX_POLY_X':                  4,
            'PRIMER_MAX_END_GC':                  2,
            'PRIMER_PRODUCT_SIZE_RANGE':          [[100, 600]],
            'PRIMER_MAX_TEMPLATE_MISPRIMING':     12.0,
            'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING':24.0,
            'PRIMER_MAX_SELF_ANY':                8.0,
            'PRIMER_MAX_SELF_END':                3.0,
            'PRIMER_PAIR_MAX_COMPL_ANY':          8.0,
            'PRIMER_PAIR_MAX_COMPL_END':          3.0,
            'PRIMER_MAX_HAIRPIN_TH':              24.0,
            'PRIMER_MAX_DIFF_TM':                 2.0,
            'PRIMER_SALT_MONOVALENT':             50.0,
            'PRIMER_SALT_DIVALENT':               1.5,
            'PRIMER_DNTP_CONC':                   0.8,
            'PRIMER_DNA_CONC':                    50.0,
        }

        self.analysis_results       = None
        self.project_name           = ""
        self.current_page           = "Analysis"
        self.input_directory        = os.getcwd()
        self.blast_family           = None
        self.blast_in_progress      = False
        self.blast_input_file       = None
        self.blast_fasta_file       = None
        self.blast_genbank_file     = None
        self.blast_gff_file         = None
        self.blast_gtf_file         = None
        self.blast_genbank_dict     = {}
        self.blast_genbank_features = {}
        self.expanded_gene_regions  = {}
        self._kmer_index: dict      = {}
        self._kmer_indices_built    = False
        self.blast_database_type    = "nt"
        self.marker_name_map        = {}
        self.reference_sequences    = {}
        self.seq_id_map             = {}
        self.id_mapping             = {}
        self.current_results_df     = None

        self.fasta_file              = tk.StringVar()
        self.excel_file              = tk.StringVar()
        self.output_file             = tk.StringVar()
        self.generate_gel_var        = tk.BooleanVar(value=True)
        self.gel_by_scaffold_var     = tk.BooleanVar(value=True)
        self.gel_by_marker_var       = tk.BooleanVar(value=True)
        self.compute_relationship_var= tk.BooleanVar(value=False)
        self.bootstrap_var           = tk.BooleanVar(value=False)
        self.bootstrap_reps          = tk.IntVar(value=100)

        self.create_primer_vars()
        self.create_widgets()

        style = ttk.Style()
        style.configure('Accent.TButton', font=('Helvetica', 10, 'bold'))

    def _poll_gui_queue(self):
        try:
            while True:
                callback = self._gui_queue.get_nowait()
                callback()
        except queue.Empty:
            pass
        self.root.after(50, self._poll_gui_queue)

    def _gui(self, callback):
        self._gui_queue.put(callback)

    def _compile_ssr_patterns(self):
        for rtype, info in self.min_repeats.items():
            motif_len  = info['motif_length']
            min_reps   = info['min_repeats']
            pattern    = r'([ATCG]{%d})\1{%d,}' % (motif_len, min_reps - 1)
            self._ssr_patterns[rtype] = re.compile(pattern, re.IGNORECASE)

    @staticmethod
    def _canonical_motif(motif: str) -> str:
        motif = motif.upper()
        rc    = ''.join(IUPAC_COMPLEMENT.get(b, b) for b in reversed(motif))
        candidates = []
        for seq in (motif, rc):
            for i in range(len(seq)):
                candidates.append(seq[i:] + seq[:i])
        return min(candidates)

    _KMER_K = 20

    def _build_kmer_index(self, seq_id: str, seq: str) -> None:
        """Constrói índice kmer para uma sequência específica"""
        k   = self._KMER_K
        idx: dict = defaultdict(list)
        for i in range(len(seq) - k + 1):
            idx[seq[i:i + k]].append(i)
        self._kmer_index[seq_id] = dict(idx)

    def _build_all_kmer_indices(self) -> None:
        """Constrói índices kmer para todas as sequências do genoma de uma vez"""
        if self._kmer_indices_built:
            return
        self.blast_status.config(text="Building kmer indices for all sequences...")
        self.root.update_idletasks()

        total = len(self.blast_genbank_dict)
        for idx, (gb_id, gb_seq) in enumerate(self.blast_genbank_dict.items()):
            if gb_id not in self._kmer_index:
                self._build_kmer_index(gb_id, gb_seq)
            if idx % 10 == 0:
                self.blast_status.config(
                    text=f"Building kmer indices: {idx+1}/{total} sequences...")
                self.root.update_idletasks()

        self._kmer_indices_built = True
        self.blast_status.config(text=f"Kmer indices built for {total} sequences")
        self.root.update_idletasks()
        logging.info("Kmer indices built for %d sequences", total)

    def _find_with_mismatch(
        self,
        region: str,
        gb_id: str,
        gb_seq: str,
        max_mismatch_rate: float,
    ) -> list:
        k        = self._KMER_K
        min_id   = (1.0 - max_mismatch_rate) * 100.0
        region_u = region.upper()
        rlen     = len(region_u)

        pos = gb_seq.upper().find(region_u)
        if pos >= 0:
            return [(pos + 1, 100.0)]

        if rlen < k:
            return []

        if gb_id not in self._kmer_index:
            self._build_kmer_index(gb_id, gb_seq)

        anchor       = region_u[:k]
        kmer_idx     = self._kmer_index[gb_id]
        candidate_starts: set = set()

        for kmer, positions in kmer_idx.items():
            mismatches = sum(a != b for a, b in zip(anchor, kmer))
            if mismatches <= max(2, int(k * max_mismatch_rate)):
                for p in positions:
                    candidate_starts.add(p)

        gb_upper = gb_seq.upper()
        results  = []
        for start in candidate_starts:
            end = start + rlen
            if end > len(gb_upper):
                continue
            segment  = gb_upper[start:end]
            matches  = sum(a == b for a, b in zip(region_u, segment))
            identity = (matches / rlen) * 100.0
            if identity >= min_id:
                results.append((start + 1, identity))

        results.sort(key=lambda x: -x[1])
        return results

    @staticmethod
    def reverse_complement(seq: str) -> str:
        return ''.join(
            IUPAC_COMPLEMENT.get(base, base) for base in reversed(seq)
        )

    def create_primer_vars(self):
        fields = [
            ('PRIMER_MIN_SIZE',          str(self.primer_params['PRIMER_MIN_SIZE'])),
            ('PRIMER_MAX_SIZE',          str(self.primer_params['PRIMER_MAX_SIZE'])),
            ('PRIMER_OPT_SIZE',          str(self.primer_params['PRIMER_OPT_SIZE'])),
            ('PRIMER_MIN_GC',            str(self.primer_params['PRIMER_MIN_GC'])),
            ('PRIMER_MAX_GC',            str(self.primer_params['PRIMER_MAX_GC'])),
            ('PRIMER_OPT_GC_PERCENT',    str(self.primer_params['PRIMER_OPT_GC_PERCENT'])),
            ('PRIMER_MIN_TM',            str(self.primer_params['PRIMER_MIN_TM'])),
            ('PRIMER_MAX_TM',            str(self.primer_params['PRIMER_MAX_TM'])),
            ('PRIMER_OPT_TM',            str(self.primer_params['PRIMER_OPT_TM'])),
            ('PRIMER_MAX_DIFF_TM',       str(self.primer_params['PRIMER_MAX_DIFF_TM'])),
            ('PRIMER_PRODUCT_SIZE_MIN',  str(self.primer_params['PRIMER_PRODUCT_SIZE_RANGE'][0][0])),
            ('PRIMER_PRODUCT_SIZE_MAX',  str(self.primer_params['PRIMER_PRODUCT_SIZE_RANGE'][0][1])),
            ('PRIMER_MAX_POLY_X',        str(self.primer_params['PRIMER_MAX_POLY_X'])),
            ('PRIMER_MAX_END_GC',        str(self.primer_params['PRIMER_MAX_END_GC'])),
            ('PRIMER_SALT_MONOVALENT',   str(self.primer_params['PRIMER_SALT_MONOVALENT'])),
            ('PRIMER_SALT_DIVALENT',     str(self.primer_params['PRIMER_SALT_DIVALENT'])),
            ('PRIMER_DNTP_CONC',         str(self.primer_params['PRIMER_DNTP_CONC'])),
            ('PRIMER_DNA_CONC',          str(self.primer_params['PRIMER_DNA_CONC'])),
            ('PRIMER_MAX_SELF_ANY',      str(self.primer_params['PRIMER_MAX_SELF_ANY'])),
            ('PRIMER_MAX_SELF_END',      str(self.primer_params['PRIMER_MAX_SELF_END'])),
            ('PRIMER_PAIR_MAX_COMPL_ANY',str(self.primer_params['PRIMER_PAIR_MAX_COMPL_ANY'])),
            ('PRIMER_PAIR_MAX_COMPL_END',str(self.primer_params['PRIMER_PAIR_MAX_COMPL_END'])),
        ]
        for name, default in fields:
            var = tk.StringVar(value=default)
            setattr(self, f"{name}_var", var)

        self._primer_entries: dict = {}

    def _attach_validation_trace(self, key: str, var: tk.StringVar,
                                 is_int: bool = False):
        def _validate(*_):
            widget = self._primer_entries.get(key)
            if widget is None:
                return
            try:
                v = int(var.get()) if is_int else float(var.get())
                if v <= 0:
                    raise ValueError
                widget.configure(foreground='black')
            except ValueError:
                widget.configure(foreground='red')

        var.trace_add('write', _validate)

    def create_widgets(self):
        header_frame = tk.Frame(self.root, bg='#2c3e50', height=100)
        header_frame.pack(fill=tk.X, pady=(0, 10))
        header_frame.pack_propagate(False)

        tk.Label(header_frame, text="SSRdev",
                 font=('Helvetica', 28, 'bold'),
                 fg='#ecf0f1', bg='#2c3e50', pady=20).pack(expand=True)
        tk.Label(header_frame,
                 text="a platform for SSR markers development — v4.5",
                 font=('Helvetica', 12), fg='#bdc3c7', bg='#2c3e50').pack()

        self.main_container = ttk.Frame(self.root)
        self.main_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.create_sidebar()
        self.create_main_content()

    def create_sidebar(self):
        sidebar = ttk.Frame(self.main_container, width=250)
        sidebar.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        sidebar.pack_propagate(False)

        scrollbar = ttk.Scrollbar(sidebar, orient="vertical")
        scrollbar.pack(side="right", fill="y")

        sb_canvas = tk.Canvas(sidebar, highlightthickness=0,
                              yscrollcommand=scrollbar.set)
        sb_canvas.pack(side="left", fill="both", expand=True)
        scrollbar.config(command=sb_canvas.yview)

        scroll_frame = ttk.Frame(sb_canvas)
        _win_id = sb_canvas.create_window((0, 0), window=scroll_frame,
                                          anchor="nw")

        def _on_canvas_resize(event):
            sb_canvas.itemconfig(_win_id, width=event.width)
        sb_canvas.bind("<Configure>", _on_canvas_resize)

        def _on_frame_configure(event):
            sb_canvas.configure(scrollregion=sb_canvas.bbox("all"))
        scroll_frame.bind("<Configure>", _on_frame_configure)

        def _on_mousewheel(event):
            # Compatibilidade entre diferentes sistemas operacionais
            if hasattr(event, 'delta'):
                sb_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
            elif event.num == 4:
                sb_canvas.yview_scroll(-1, "units")
            elif event.num == 5:
                sb_canvas.yview_scroll(1, "units")

        def _bind_scroll(widget):
            widget.bind("<MouseWheel>", _on_mousewheel, add="+")
            widget.bind("<Button-4>",   _on_mousewheel, add="+")
            widget.bind("<Button-5>",   _on_mousewheel, add="+")
            for child in widget.winfo_children():
                _bind_scroll(child)

        sb_canvas.bind("<MouseWheel>", _on_mousewheel)
        sb_canvas.bind("<Button-4>",   _on_mousewheel)
        sb_canvas.bind("<Button-5>",   _on_mousewheel)

        logo_frame = ttk.Frame(scroll_frame)
        logo_frame.pack(fill=tk.X, padx=10, pady=10)
        ttk.Label(logo_frame, text="🧬",
                  font=('Helvetica', 24)).pack()
        ttk.Label(logo_frame, text="SSR Marker Development",
                  font=('Helvetica', 10, 'bold')).pack()

        ttk.Separator(scroll_frame, orient=tk.HORIZONTAL).pack(
            fill=tk.X, padx=10, pady=5)

        nav_frame = ttk.Frame(scroll_frame)
        nav_frame.pack(fill=tk.X, padx=10, pady=10)
        ttk.Label(nav_frame, text="NAVIGATION",
                  font=('Helvetica', 9, 'bold'),
                  foreground='#7f8c8d').pack(anchor=tk.W)

        for option in ["Prospecting", "Primer Settings",
                       "Characterization", "Validation"]:
            ttk.Button(nav_frame, text=option,
                       command=lambda o=option: self.show_page(o),
                       style="Nav.TButton").pack(fill=tk.X, pady=2)

        ttk.Separator(scroll_frame, orient=tk.HORIZONTAL).pack(
            fill=tk.X, padx=10, pady=5)
        proj_frame = ttk.LabelFrame(scroll_frame, text="Project", padding=8)
        proj_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Button(proj_frame, text="💾  Save Project",
                   command=self.save_project).pack(fill=tk.X, pady=2)
        ttk.Button(proj_frame, text="📂  Load Project",
                   command=self.load_project).pack(fill=tk.X, pady=2)

        ttk.Separator(scroll_frame, orient=tk.HORIZONTAL).pack(
            fill=tk.X, padx=10, pady=5)

        config_frame = ttk.LabelFrame(scroll_frame,
                                      text="SSR Search Settings", padding=10)
        config_frame.pack(fill=tk.X, padx=10, pady=5)

        for rtype, info in self.min_repeats.items():
            row = ttk.Frame(config_frame)
            row.pack(fill=tk.X, pady=3)

            enabled_var = tk.BooleanVar(value=info['enabled'])
            ttk.Checkbutton(row, variable=enabled_var).pack(side=tk.LEFT)
            setattr(self, f"{rtype}_enabled_var", enabled_var)

            ttk.Label(row,
                      text=f"{rtype.capitalize()} ({info['motif_length']}bp):"
                      ).pack(side=tk.LEFT, padx=(2, 4))

            var = tk.StringVar(value=str(info['min_repeats']))
            ttk.Spinbox(row, from_=1, to=100,
                        textvariable=var, width=4).pack(side=tk.RIGHT)
            setattr(self, f"{rtype}_var", var)

        ttk.Separator(scroll_frame, orient=tk.HORIZONTAL).pack(
            fill=tk.X, padx=10, pady=8)
        ttk.Label(scroll_frame, text="Version 4.5",
                  font=('Helvetica', 8),
                  foreground='#7f8c8d').pack(pady=(0, 10))

        self.root.after(100, lambda: _bind_scroll(scroll_frame))

        self.root.style = ttk.Style()
        self.root.style.configure('Nav.TButton',
                                  anchor=tk.W, padding=(10, 5))

    def create_main_content(self):
        self.notebook = ttk.Notebook(self.main_container)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        self.analysis_page = ttk.Frame(self.notebook)
        self.setup_analysis_page()
        self.notebook.add(self.analysis_page, text="Prospecting")

        self.primer_page = ttk.Frame(self.notebook)
        self.setup_primer_page()
        self.notebook.add(self.primer_page, text="Primer Settings")

        self.characterization_page = ttk.Frame(self.notebook)
        self.setup_characterization_page()
        self.notebook.add(self.characterization_page,
                          text="Characterization")

        self.validation_page = ttk.Frame(self.notebook)
        self.setup_validation_page()
        self.notebook.add(self.validation_page, text="Validation")

    def setup_analysis_page(self):
        main_frame = ttk.Frame(self.analysis_page, padding=20)
        main_frame.pack(fill=tk.BOTH, expand=True)

        proj_frame = ttk.Frame(main_frame)
        proj_frame.pack(fill=tk.X, pady=(0, 15))
        ttk.Label(proj_frame, text="Project Name:",
                  font=('Helvetica', 10, 'bold')).pack(anchor=tk.W)
        self.project_entry = ttk.Entry(proj_frame, font=('Helvetica', 10))
        self.project_entry.pack(fill=tk.X, pady=5)

        file_frame = ttk.LabelFrame(main_frame, text="FASTA Input", padding=10)
        file_frame.pack(fill=tk.X, pady=(0, 15))
        ttk.Label(file_frame, text="FASTA File:").grid(
            row=0, column=0, sticky=tk.W, pady=5)
        self.file_path = tk.StringVar()
        ttk.Entry(file_frame, textvariable=self.file_path,
                  state='readonly', width=50).grid(
            row=0, column=1, sticky=tk.EW, padx=5)
        ttk.Button(file_frame, text="Browse",
                   command=self.browse_file).grid(row=0, column=2, padx=5)
        file_frame.columnconfigure(1, weight=1)

        paste_frame = ttk.LabelFrame(main_frame,
                                     text="Or paste sequence", padding=10)
        paste_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))
        self.seq_text = scrolledtext.ScrolledText(
            paste_frame, height=10, font=('Courier', 10))
        self.seq_text.pack(fill=tk.BOTH, expand=True)

        progress_frame = ttk.Frame(main_frame)
        progress_frame.pack(fill=tk.X, pady=(0, 15))
        self.progress = ttk.Progressbar(
            progress_frame, orient=tk.HORIZONTAL, mode='determinate')
        self.progress.pack(fill=tk.X)
        self.status = ttk.Label(progress_frame, text="Ready",
                                background='#f0f0f0')
        self.status.pack(fill=tk.X, pady=5)

        btn_frame = ttk.Frame(main_frame)
        btn_frame.pack(fill=tk.X, pady=(0, 15))
        button_container = ttk.Frame(btn_frame)
        button_container.pack(expand=True, pady=5)

        self.run_analysis_btn = ttk.Button(
            button_container, text="Run Analysis",
            command=self.run_analysis, style="Accent.TButton")
        self.run_analysis_btn.pack(side=tk.LEFT, padx=5)

        ttk.Button(button_container, text="Clear Results",
                   command=self.clear_results).pack(side=tk.LEFT, padx=5)

        results_outer = ttk.LabelFrame(main_frame, text="Results", padding=10)
        results_outer.pack(fill=tk.BOTH, expand=True)
        self.results_frame = ttk.Frame(results_outer)
        self.results_frame.pack(fill=tk.BOTH, expand=True)

    def setup_primer_page(self):
        main_canvas   = tk.Canvas(self.primer_page, borderwidth=0)
        scrollbar     = ttk.Scrollbar(self.primer_page, orient=tk.VERTICAL,
                                      command=main_canvas.yview)
        scroll_frame  = ttk.Frame(main_canvas)

        scroll_frame.bind(
            "<Configure>",
            lambda e: main_canvas.configure(
                scrollregion=main_canvas.bbox("all"))
        )
        main_canvas.create_window((0, 0), window=scroll_frame, anchor="nw")
        main_canvas.configure(yscrollcommand=scrollbar.set)
        main_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        params = [
            ('Primer Size (bp)',           'PRIMER_MIN_SIZE', 'PRIMER_MAX_SIZE', 'PRIMER_OPT_SIZE',          True),
            ('GC Content (%)',             'PRIMER_MIN_GC',   'PRIMER_MAX_GC',   'PRIMER_OPT_GC_PERCENT',    False),
            ('Melting Temp (°C)',          'PRIMER_MIN_TM',   'PRIMER_MAX_TM',   'PRIMER_OPT_TM',            False),
            ('Max Tm Difference (°C)',     'PRIMER_MAX_DIFF_TM', None, None,                                  False),
            ('Product Size Range (bp)',    'PRIMER_PRODUCT_SIZE_MIN', 'PRIMER_PRODUCT_SIZE_MAX', None,         True),
            ('Max Poly-X Length',          'PRIMER_MAX_POLY_X', None, None,                                   True),
            ('Max 3\' GC Clamp',           'PRIMER_MAX_END_GC', None, None,                                   True),
            ('Salt Concentration (mM)',    'PRIMER_SALT_MONOVALENT', None, None,                              False),
            ('Divalent Concentration (mM)','PRIMER_SALT_DIVALENT', None, None,                                False),
            ('dNTP Concentration (mM)',    'PRIMER_DNTP_CONC', None, None,                                    False),
            ('DNA Concentration (nM)',     'PRIMER_DNA_CONC', None, None,                                     False),
            ('Max Self Complementarity',   'PRIMER_MAX_SELF_ANY', None, None,                                 False),
            ('Max 3\' Self Complementarity','PRIMER_MAX_SELF_END', None, None,                                False),
            ('Max Pair Complementarity',   'PRIMER_PAIR_MAX_COMPL_ANY', None, None,                           False),
            ('Max 3\' Pair Complementarity','PRIMER_PAIR_MAX_COMPL_END', None, None,                          False),
        ]

        for label, min_key, max_key, opt_key, is_int in params:
            frame = ttk.LabelFrame(scroll_frame, text=label, padding=10)
            frame.pack(fill=tk.X, padx=10, pady=5)
            col = 0
            for role, key in [('Min:', min_key), ('Opt:', opt_key),
                               ('Max:', max_key)]:
                if key is None:
                    continue
                ttk.Label(frame, text=role).grid(row=0, column=col, padx=5)
                var    = getattr(self, f"{key}_var")
                entry  = ttk.Entry(frame, textvariable=var, width=8)
                entry.grid(row=0, column=col + 1, padx=5)
                self._primer_entries[key] = entry
                self._attach_validation_trace(key, var, is_int)
                col += 2

        btn_frame = ttk.Frame(scroll_frame)
        btn_frame.pack(fill=tk.X, padx=10, pady=10)
        ttk.Button(btn_frame, text="Save Settings",
                   command=self.save_primer_settings,
                   style="Accent.TButton").pack(side=tk.LEFT)
        ttk.Button(btn_frame, text="Reset to Defaults",
                   command=self.reset_primer_settings).pack(
            side=tk.LEFT, padx=5)

    def setup_characterization_page(self):
        main_canvas  = tk.Canvas(self.characterization_page, borderwidth=0)
        scrollbar    = ttk.Scrollbar(self.characterization_page,
                                     orient=tk.VERTICAL,
                                     command=main_canvas.yview)
        scroll_frame = ttk.Frame(main_canvas)

        scroll_frame.bind(
            "<Configure>",
            lambda e: main_canvas.configure(
                scrollregion=main_canvas.bbox("all"))
        )
        main_canvas.create_window((0, 0), window=scroll_frame, anchor="nw")
        main_canvas.configure(yscrollcommand=scrollbar.set)
        main_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        main_frame = ttk.Frame(scroll_frame, padding=20)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # ========== FAMILY SELECTION ==========
        family_frame = ttk.LabelFrame(
            main_frame, text="NCBI BLAST Family Selection", padding=10)
        family_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(family_frame, text="Select family:").grid(
            row=0, column=0, sticky=tk.W, pady=5)
        self.family_var = tk.StringVar()
        ttk.Entry(family_frame, textvariable=self.family_var,
                  state='readonly', width=40).grid(
            row=0, column=1, sticky=tk.EW, padx=5)
        ttk.Button(family_frame, text="Choose Family",
                   command=self.select_family).grid(
            row=0, column=2, padx=5)
        family_frame.columnconfigure(1, weight=1)

        ttk.Label(family_frame, text="BLAST Database:").grid(
            row=1, column=0, sticky=tk.W, pady=5)
        self.blast_db_var = tk.StringVar(value="DNA_refSeq")
        db_frame = ttk.Frame(family_frame)
        db_frame.grid(row=1, column=1, sticky=tk.W, padx=5)
        ttk.Radiobutton(db_frame, text="DNA_refSeq", variable=self.blast_db_var,
                        value="DNA_refSeq").pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(db_frame, text="mRNA_refSeq", variable=self.blast_db_var,
                        value="mRNA_refSeq").pack(side=tk.LEFT, padx=5)

        ttk.Separator(family_frame, orient=tk.HORIZONTAL).grid(
            row=2, column=0, columnspan=3, sticky=tk.EW, pady=10)

        # ========== LOCAL GENBANK ==========
        ttk.Label(family_frame,
                  text="OR use Local GenBank Reference Genome:",
                  font=('Helvetica', 10, 'bold')).grid(
            row=3, column=0, columnspan=3, sticky=tk.W, pady=5)
        ttk.Label(family_frame, text="GenBank File:").grid(
            row=4, column=0, sticky=tk.W, pady=5)
        self.genbank_var = tk.StringVar()
        ttk.Entry(family_frame, textvariable=self.genbank_var,
                  state='readonly', width=40).grid(
            row=4, column=1, sticky=tk.EW, padx=5)
        ttk.Button(family_frame, text="Browse GenBank",
                   command=self.browse_genbank).grid(row=4, column=2, padx=5)
        ttk.Label(family_frame,
                  text="GenBank files can be downloaded from NCBI "
                       "(format: .gb, .gbff, .gbk)").grid(
            row=5, column=0, columnspan=3, sticky=tk.W, pady=5)

        # ========== GFF3 ==========
        ttk.Label(family_frame,
                  text="OR use GFF3 + FASTA:",
                  font=('Helvetica', 10, 'bold')).grid(
            row=6, column=0, columnspan=3, sticky=tk.W, pady=5)
        ttk.Label(family_frame, text="GFF3 File:").grid(
            row=7, column=0, sticky=tk.W, pady=5)
        self.gff_var = tk.StringVar()
        ttk.Entry(family_frame, textvariable=self.gff_var,
                  state='readonly', width=40).grid(
            row=7, column=1, sticky=tk.EW, padx=5)
        ttk.Button(family_frame, text="Browse GFF3",
                   command=self.browse_gff).grid(row=7, column=2, padx=5)

        ttk.Label(family_frame, text="GFF3 FASTA File:").grid(
            row=8, column=0, sticky=tk.W, pady=5)
        self.gff_fasta_var = tk.StringVar()
        ttk.Entry(family_frame, textvariable=self.gff_fasta_var,
                  state='readonly', width=40).grid(
            row=8, column=1, sticky=tk.EW, padx=5)
        ttk.Button(family_frame, text="Browse FASTA for GFF3",
                   command=self.browse_gff_fasta).grid(row=8, column=2, padx=5)

        # ========== GTF ==========
        ttk.Label(family_frame,
                  text="OR use GTF + FASTA:",
                  font=('Helvetica', 10, 'bold')).grid(
            row=9, column=0, columnspan=3, sticky=tk.W, pady=5)
        ttk.Label(family_frame, text="GTF File:").grid(
            row=10, column=0, sticky=tk.W, pady=5)
        self.gtf_var = tk.StringVar()
        ttk.Entry(family_frame, textvariable=self.gtf_var,
                  state='readonly', width=40).grid(
            row=10, column=1, sticky=tk.EW, padx=5)
        ttk.Button(family_frame, text="Browse GTF",
                   command=self.browse_gtf).grid(row=10, column=2, padx=5)

        ttk.Label(family_frame, text="GTF FASTA File:").grid(
            row=11, column=0, sticky=tk.W, pady=5)
        self.gtf_fasta_var = tk.StringVar()
        ttk.Entry(family_frame, textvariable=self.gtf_fasta_var,
                  state='readonly', width=40).grid(
            row=11, column=1, sticky=tk.EW, padx=5)
        ttk.Button(family_frame, text="Browse FASTA for GTF",
                   command=self.browse_gtf_fasta).grid(row=11, column=2, padx=5)

        btn_frame_annot = ttk.Frame(family_frame)
        btn_frame_annot.grid(row=12, column=0, columnspan=3, pady=10)
        ttk.Button(btn_frame_annot, text="Verify Annotations",
                   command=self.verify_genbank_annotations).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame_annot, text="Clear GenBank Data",
                   command=self.clear_genbank_data).pack(side=tk.LEFT, padx=5)

        # ========== REFERENCE FASTA ==========
        fasta_ref_frame = ttk.LabelFrame(
            main_frame, text="Reference FASTA File (from Prospecting)", padding=10)
        fasta_ref_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Label(fasta_ref_frame, text="Select FASTA file:").grid(
            row=0, column=0, sticky=tk.W, pady=5)
        self.blast_fasta_var = tk.StringVar()
        ttk.Entry(fasta_ref_frame, textvariable=self.blast_fasta_var,
                  state='readonly', width=40).grid(
            row=0, column=1, sticky=tk.EW, padx=5)
        ttk.Button(fasta_ref_frame, text="Browse FASTA",
                   command=self.browse_blast_fasta).grid(
            row=0, column=2, padx=5)
        ttk.Label(fasta_ref_frame,
                  text="Load the same FASTA file used during prospecting\n"
                       "This file contains the reference sequences for the SSRs",
                  font=('Helvetica', 8), foreground='gray').grid(
            row=1, column=0, columnspan=3, sticky=tk.W, pady=5)
        fasta_ref_frame.columnconfigure(1, weight=1)

        # ========== SSR RESULTS FILE ==========
        input_frame = ttk.LabelFrame(
            main_frame, text="SSR Results File", padding=10)
        input_frame.pack(fill=tk.X, pady=(0, 15))
        ttk.Label(input_frame, text="Select SSR_Results.xlsx file:").grid(
            row=0, column=0, sticky=tk.W, pady=5)
        self.blast_input_var = tk.StringVar()
        ttk.Entry(input_frame, textvariable=self.blast_input_var,
                  state='readonly', width=40).grid(
            row=0, column=1, sticky=tk.EW, padx=5)
        ttk.Button(input_frame, text="Browse",
                   command=self.browse_blast_input).grid(
            row=0, column=2, padx=5)
        input_frame.columnconfigure(1, weight=1)

        # ========== BLAST PARAMETERS ==========
        param_frame = ttk.LabelFrame(
            main_frame, text="BLAST Parameters", padding=10)
        param_frame.pack(fill=tk.X, pady=(0, 15))

        for row, (label, attr, default) in enumerate([
            ("Flanking region size (bp):", 'flank_size_var', "2000"),
            ("Gene flanking region (bp):", 'gene_flank_var', "2000"),
            ("Min overlap (bp):", 'min_overlap_var', "50"),
            ("Max mismatches (0-100%):", 'max_mismatch_var', "10"),
        ]):
            ttk.Label(param_frame, text=label).grid(
                row=row, column=0, sticky=tk.W, pady=5)
            var = tk.StringVar(value=default)
            setattr(self, attr, var)
            ttk.Entry(param_frame, textvariable=var, width=10).grid(
                row=row, column=1, padx=5)
        ttk.Label(param_frame, text="%").grid(row=3, column=2, sticky=tk.W)

        # ========== BLAST CONTROLS ==========
        controls_frame = ttk.LabelFrame(main_frame, text="BLAST Controls", padding=10)
        controls_frame.pack(fill=tk.X, pady=(0, 10))

        progress_frame = ttk.Frame(controls_frame)
        progress_frame.pack(fill=tk.X, pady=(0, 5))
        self.blast_progress = ttk.Progressbar(
            progress_frame, orient=tk.HORIZONTAL, mode='determinate')
        self.blast_progress.pack(fill=tk.X)

        self.blast_status = ttk.Label(controls_frame, text="Ready for BLAST analysis")
        self.blast_status.pack(fill=tk.X, pady=5)

        btn_frame = ttk.Frame(controls_frame)
        btn_frame.pack(fill=tk.X, pady=5)
        button_container = ttk.Frame(btn_frame)
        button_container.pack(expand=True)

        self.run_ncbi_blast_btn = ttk.Button(
            button_container, text="Run NCBI BLAST",
            command=self.run_blast_analysis,
            style="Accent.TButton")
        self.run_ncbi_blast_btn.pack(side=tk.LEFT, padx=5)

        self.run_local_blast_btn = ttk.Button(
            button_container, text="Run Local GenBank BLAST",
            command=self.run_local_blast_analysis,
            style="Accent.TButton")
        self.run_local_blast_btn.pack(side=tk.LEFT, padx=5)

        self.cancel_blast_btn = ttk.Button(
            button_container, text="Cancel BLAST",
            command=self.cancel_blast)
        self.cancel_blast_btn.pack(side=tk.LEFT, padx=5)

        # ========== BLAST RESULTS ==========
        results_container = ttk.LabelFrame(main_frame, text="BLAST Results", padding=10)
        results_container.pack(fill=tk.BOTH, expand=True)

        # Frame que será limpo
        self.blast_results_frame = ttk.Frame(results_container)
        self.blast_results_frame.pack(fill=tk.BOTH, expand=True)

        # Mensagem inicial
        ttk.Label(self.blast_results_frame,
                  text="Run BLAST to see results here",
                  font=('Helvetica', 10, 'italic'),
                  foreground='gray').pack(expand=True)

    def setup_validation_page(self):
        main_canvas  = tk.Canvas(self.validation_page, borderwidth=0)
        scrollbar    = ttk.Scrollbar(self.validation_page, orient=tk.VERTICAL,
                                     command=main_canvas.yview)
        scroll_frame = ttk.Frame(main_canvas)

        scroll_frame.bind(
            "<Configure>",
            lambda e: main_canvas.configure(
                scrollregion=main_canvas.bbox("all"))
        )
        main_canvas.create_window((0, 0), window=scroll_frame, anchor="nw")
        main_canvas.configure(yscrollcommand=scrollbar.set)
        main_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        main_frame = ttk.Frame(scroll_frame, padding=20)
        main_frame.pack(fill=tk.BOTH, expand=True)

        input_frame = ttk.LabelFrame(main_frame, text="Input Files", padding=10)
        input_frame.pack(fill=tk.X, pady=(0, 15))

        for row, (label, var_attr, browse_cmd) in enumerate([
            ("Multi-FASTA File:",   'fasta_file',  'browse_fasta'),
            ("Primer Excel File:",  'excel_file',  'browse_excel'),
        ]):
            ttk.Label(input_frame, text=label).grid(
                row=row, column=0, sticky=tk.W, pady=5)
            ttk.Entry(input_frame,
                      textvariable=getattr(self, var_attr),
                      width=50).grid(row=row, column=1,
                                     sticky=tk.EW, padx=5)
            ttk.Button(input_frame, text="Browse",
                       command=getattr(self, browse_cmd)).grid(
                row=row, column=2, padx=5)
        input_frame.columnconfigure(1, weight=1)

        options_frame = ttk.LabelFrame(main_frame, text="Options", padding=10)
        options_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Checkbutton(
            options_frame,
            text="Generate virtual gel electrophoresis images",
            variable=self.generate_gel_var).pack(anchor=tk.W, pady=2)

        gel_options = ttk.Frame(options_frame)
        gel_options.pack(fill=tk.X, padx=20, pady=5)
        ttk.Checkbutton(
            gel_options,
            text="Gel for each scaffold (all markers per scaffold)",
            variable=self.gel_by_scaffold_var).pack(anchor=tk.W, pady=2)
        ttk.Checkbutton(
            gel_options,
            text="Gel for each marker (all scaffolds per marker)",
            variable=self.gel_by_marker_var).pack(anchor=tk.W, pady=2)

        ttk.Checkbutton(
            options_frame, text="Compute relationship",
            variable=self.compute_relationship_var,
            command=self.toggle_relationship_options).pack(
            anchor=tk.W, pady=2)

        self.relationship_frame = ttk.Frame(options_frame)
        ttk.Checkbutton(self.relationship_frame,
                        text="Bootstrap dendrogram",
                        variable=self.bootstrap_var).pack(side=tk.LEFT)
        ttk.Label(self.relationship_frame,
                  text="Replicates:").pack(side=tk.LEFT, padx=(10, 5))
        ttk.Entry(self.relationship_frame,
                  textvariable=self.bootstrap_reps,
                  width=5).pack(side=tk.LEFT)
        self.relationship_frame.pack_forget()

        progress_frame = ttk.Frame(main_frame)
        progress_frame.pack(fill=tk.X, pady=(0, 15))
        self.validation_progress = ttk.Progressbar(
            progress_frame, orient=tk.HORIZONTAL, mode='determinate')
        self.validation_progress.pack(fill=tk.X)
        self.validation_status = ttk.Label(progress_frame, text="Ready")
        self.validation_status.pack(fill=tk.X, pady=5)

        btn_frame = ttk.Frame(main_frame)
        btn_frame.pack(fill=tk.X, pady=(0, 15))
        button_container = ttk.Frame(btn_frame)
        button_container.pack(expand=True)
        ttk.Button(button_container, text="Analyze Primers",
                   command=self.validate_primers,
                   style="Accent.TButton").pack(pady=10)

    # ==================== MÉTODOS DE BROWSE ====================

    def browse_blast_fasta(self):
        fp = filedialog.askopenfilename(
            filetypes=[('FASTA files', '*.fasta *.fa *.fna'),
                       ('All files', '*.*')])
        if fp:
            self.blast_fasta_var.set(fp)
            self.blast_fasta_file = fp
            self._load_reference_sequences(fp)

    def _load_reference_sequences(self, fasta_file: str):
        try:
            self.blast_status.config(
                text=f"Loading reference sequences from: {os.path.basename(fasta_file)}…")
            self.root.update_idletasks()

            self.reference_sequences = {}
            self.seq_id_map = {}
            for record in SeqIO.parse(fasta_file, "fasta"):
                original_id = record.id
                seq_str = str(record.seq).upper()
                self.reference_sequences[original_id] = seq_str

                self.seq_id_map[original_id] = original_id

                if '_' in original_id:
                    first_part = original_id.split('_')[0]
                    if first_part not in self.seq_id_map:
                        self.seq_id_map[first_part] = original_id

                clean_id = re.sub(r'[^\w]', '_', original_id)
                if clean_id != original_id and clean_id not in self.seq_id_map:
                    self.seq_id_map[clean_id] = original_id

            self.blast_status.config(
                text=f"Loaded {len(self.reference_sequences)} reference sequences")
            logging.info("Reference sequences loaded: %s sequences from %s",
                        len(self.reference_sequences), fasta_file)
        except Exception as e:
            error_msg = str(e)
            self.blast_status.config(text="Error loading reference sequences")
            messagebox.showerror("Error", f"Failed to load FASTA file:\n{error_msg}")
            logging.error("FASTA load error: %s", error_msg)

    def browse_gff(self):
        fp = filedialog.askopenfilename(
            filetypes=[('GFF3 files', '*.gff *.gff3'), ('All files', '*.*')])
        if fp:
            self.gff_var.set(fp)
            self.blast_gff_file = fp

    def browse_gff_fasta(self):
        fp = filedialog.askopenfilename(
            filetypes=[('FASTA files', '*.fasta *.fa *.fna'), ('All files', '*.*')])
        if fp:
            self.gff_fasta_var.set(fp)
            self.load_gff_annotations(self.blast_gff_file, fp)

    def browse_gtf(self):
        fp = filedialog.askopenfilename(
            filetypes=[('GTF files', '*.gtf'), ('All files', '*.*')])
        if fp:
            self.gtf_var.set(fp)
            self.blast_gtf_file = fp

    def browse_gtf_fasta(self):
        fp = filedialog.askopenfilename(
            filetypes=[('FASTA files', '*.fasta *.fa *.fna'), ('All files', '*.*')])
        if fp:
            self.gtf_fasta_var.set(fp)
            self.load_gtf_annotations(self.blast_gtf_file, fp)

    def _load_annotation_common(self, annotation_file: str, fasta_file: str,
                                annotation_type: str, parse_func):
        try:
            self.blast_status.config(
                text=f"Loading {annotation_type} annotations from {os.path.basename(annotation_file)}…")
            self.root.update_idletasks()

            self.blast_genbank_dict = {}
            self.blast_genbank_features = {}
            self._kmer_indices_built = False

            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_id = record.id
                seq_str = str(record.seq).upper()
                self.blast_genbank_dict[seq_id] = seq_str
                self.blast_genbank_features[seq_id] = []

            features = parse_func(annotation_file)

            for feature in features:
                seq_id = feature['seqid']
                if seq_id not in self.blast_genbank_features:
                    self.blast_genbank_features[seq_id] = []

                self.blast_genbank_features[seq_id].append({
                    'type': feature['type'],
                    'start': feature['start'],
                    'end': feature['end'],
                    'gene': feature['gene'],
                    'strand': feature['strand']
                })

            n = len(self.blast_genbank_dict)
            total_genes = sum(len(f) for f in self.blast_genbank_features.values())
            self.blast_status.config(
                text=f"Loaded {n} sequences with {total_genes} annotated features from {annotation_type}")
            messagebox.showinfo("Success",
                f"{annotation_type} file loaded!\nFound {n} sequences with {total_genes} annotated features.")
            logging.info("%s loaded: %s sequences from %s", annotation_type, n, annotation_file)

        except Exception as e:
            error_msg = str(e)
            messagebox.showerror("Error", f"Failed to load {annotation_type} file:\n{error_msg}")
            self.blast_status.config(text=f"Error loading {annotation_type} file")
            logging.error("%s load error: %s", annotation_type, error_msg)

    def load_gff_annotations(self, gff_file: str, fasta_file: str):
        def parse_gff(gff_file):
            db = gffutils.create_db(gff_file, dbfn=':memory:', force=True,
                                    keep_order=True, merge_strategy='merge')
            features = []
            for feature in db.features_of_type(['gene', 'CDS', 'mRNA', 'tRNA', 'rRNA']):
                gene_name = feature.attributes.get('Name',
                            feature.attributes.get('ID', ['unknown'])[0])
                features.append({
                    'seqid': feature.seqid,
                    'type': feature.featuretype,
                    'start': feature.start,
                    'end': feature.end,
                    'gene': str(gene_name),
                    'strand': feature.strand if feature.strand else 1
                })
            return features

        self._load_annotation_common(gff_file, fasta_file, "GFF3", parse_gff)

    def load_gtf_annotations(self, gtf_file: str, fasta_file: str):
        def parse_gtf(gtf_file):
            df = gtfparse.read_gtf(gtf_file)
            features = []
            for _, row in df.iterrows():
                if row['feature'] in ['gene', 'transcript', 'exon', 'CDS']:
                    gene_name = row.get('gene_name', row.get('gene_id', 'unknown'))
                    features.append({
                        'seqid': row['seqname'],
                        'type': row['feature'],
                        'start': row['start'],
                        'end': row['end'],
                        'gene': str(gene_name),
                        'strand': 1 if row['strand'] == '+' else -1
                    })
            return features

        self._load_annotation_common(gtf_file, fasta_file, "GTF", parse_gtf)

    def verify_genbank_annotations(self):
        if not self.blast_genbank_features and not self.blast_genbank_dict:
            if self.blast_genbank_file:
                self.browse_genbank()
            elif self.blast_gff_file:
                if self.gff_fasta_var.get():
                    self.load_gff_annotations(self.blast_gff_file, self.gff_fasta_var.get())
                else:
                    messagebox.showwarning("Warning", "Please select both GFF3 and FASTA files first")
                    return
            elif self.blast_gtf_file:
                if self.gtf_fasta_var.get():
                    self.load_gtf_annotations(self.blast_gtf_file, self.gtf_fasta_var.get())
                else:
                    messagebox.showwarning("Warning", "Please select both GTF and FASTA files first")
                    return
            else:
                messagebox.showwarning("Warning", "No annotation file loaded")
                return

        stats = {}
        for seq_id, features in self.blast_genbank_features.items():
            stats[seq_id] = len(features)

        total_genes = sum(stats.values())
        sequences_with_annotations = sum(1 for v in stats.values() if v > 0)

        messagebox.showinfo("Annotations Summary",
            f"Loaded {len(self.blast_genbank_dict)} sequences\n"
            f"Total annotated features: {total_genes}\n"
            f"Sequences with annotations: {sequences_with_annotations}\n"
            f"Sequences without annotations: {len(self.blast_genbank_dict) - sequences_with_annotations}\n\n"
            f"First 5 sequences with annotations:\n" +
            "\n".join([f"  {seq_id}: {count} features" for seq_id, count in list(stats.items())[:5] if count > 0]))

    def clear_genbank_data(self):
        self.blast_genbank_dict = {}
        self.blast_genbank_features = {}
        self.expanded_gene_regions = {}
        self._kmer_index = {}
        self._kmer_indices_built = False
        self.id_mapping = {}
        self.blast_genbank_file = None
        self.blast_gff_file = None
        self.blast_gtf_file = None
        self.genbank_var.set("")
        self.gff_var.set("")
        self.gtf_var.set("")
        self.gff_fasta_var.set("")
        self.gtf_fasta_var.set("")
        self.blast_status.config(text="GenBank/GFF/GTF data cleared")
        messagebox.showinfo("Info", "GenBank/GFF/GTF data has been cleared.\nYou can now load a new annotation file.")

    def browse_genbank(self):
        fp = filedialog.askopenfilename(
            filetypes=[('GenBank files', '*.gb *.gbff *.gbk'),
                       ('All files', '*.*')])
        if not fp:
            return

        self.genbank_var.set(fp)
        self.blast_genbank_file = fp

        try:
            self.blast_status.config(
                text=f"Loading GenBank file: {os.path.basename(fp)}…")
            self.root.update_idletasks()

            self.blast_genbank_dict     = {}
            self.blast_genbank_features = {}
            self.expanded_gene_regions  = {}
            self._kmer_index            = {}
            self._kmer_indices_built    = False

            records = list(SeqIO.parse(fp, "genbank"))
            self.id_mapping = self.create_id_mapping(records)

            for record in records:
                seq_id  = record.id
                seq_str = str(record.seq).upper()
                self.blast_genbank_dict[seq_id]     = seq_str
                self.blast_genbank_features[seq_id] = []
                self.expanded_gene_regions[seq_id]  = []

                if (hasattr(record, 'name') and record.name
                        and record.name != seq_id):
                    self.blast_genbank_dict[record.name] = seq_str
                    self.blast_genbank_features[record.name] = \
                        self.blast_genbank_features[seq_id]

                if hasattr(record, 'features'):
                    for feature in record.features:
                        if feature.type in ['gene', 'CDS', 'mRNA',
                                            'tRNA', 'rRNA']:
                            if (hasattr(feature.location, 'start') and
                                    hasattr(feature.location, 'end')):
                                gene_name = (
                                    feature.qualifiers.get('gene', ['unknown'])[0]
                                    if 'gene' in feature.qualifiers
                                    else feature.qualifiers.get(
                                        'product', ['unknown'])[0]
                                )
                                fi = {
                                    'type':   feature.type,
                                    'start':  int(feature.location.start) + 1,
                                    'end':    int(feature.location.end),
                                    'gene':   gene_name,
                                    'strand': getattr(
                                        feature.location, 'strand', 1),
                                }
                                self.blast_genbank_features[seq_id].append(fi)

            has_genes = any(len(f) > 0 for f in self.blast_genbank_features.values())
            if not has_genes:
                messagebox.showwarning("Warning",
                    "GenBank file has no gene annotations. "
                    "Try downloading a 'genomic.gbff' file from NCBI "
                    "which includes annotated features.")

            n = len(self.blast_genbank_dict)
            total_genes = sum(len(f) for f in self.blast_genbank_features.values())
            self.blast_status.config(
                text=f"Loaded {n} sequences with {total_genes} annotated features from GenBank")
            messagebox.showinfo(
                "Success",
                f"GenBank file loaded!\nFound {n} sequences with {total_genes} annotated features.")
            logging.info("GenBank loaded: %s sequences from %s", n, fp)

        except Exception as gb_exc:
            error_msg = str(gb_exc)
            messagebox.showerror("Error",
                                 f"Failed to load GenBank file:\n{error_msg}")
            self.blast_genbank_dict     = {}
            self.blast_genbank_features = {}
            self.blast_status.config(text="Error loading GenBank file")
            logging.error("GenBank load error: %s", error_msg)

    def create_id_mapping(self, genbank_records):
        mapping = {}
        for record in genbank_records:
            mapping[record.id] = record.id
            if hasattr(record, 'name') and record.name:
                mapping[record.name] = record.id
            if hasattr(record, 'description') and record.description:
                desc_id = record.description.split()[0]
                mapping[desc_id] = record.id
            if '.' in record.id:
                base_id = record.id.split('.')[0]
                mapping[base_id] = record.id
        return mapping

    # ==================== MÉTODOS DE ANÁLISE ====================

    def run_analysis(self):
        if self.analysis_in_progress:
            messagebox.showwarning("Busy", "An analysis is already running.")
            return

        sequences = self.get_sequences()
        if not sequences:
            return

        if not self.update_parameters():
            return

        self.project_name = self.project_entry.get().strip()

        self.analysis_in_progress = True
        self.run_analysis_btn.config(state='disabled')
        self.analysis_results = []
        self.progress['maximum'] = len(sequences)
        self.progress['value']   = 0
        self.status.config(text="Starting analysis…")

        threading.Thread(target=self._run_analysis_thread,
                         args=(sequences,), daemon=True).start()

    def _run_analysis_thread(self, sequences):
        results = []

        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = {
                executor.submit(self._analyse_one_sequence, record): record
                for record in sequences
            }
            done = 0
            for future in concurrent.futures.as_completed(futures):
                done += 1
                record = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    logging.error("Error analysing %s: %s", record.id, exc)
                    result = {'id': record.id,
                              'sequence': str(record.seq),
                              'results': []}
                results.append(result)

                _done = done
                _id   = record.id
                self._gui(lambda d=_done, i=_id: (
                    self.progress.__setitem__('value', d),
                    self.status.config(text=f"Analysed {i} ({d}/{len(sequences)})")
                ))

        order = {r.id: idx for idx, r in enumerate(sequences)}
        results.sort(key=lambda x: order.get(x['id'], 999999))
        self.analysis_results = results

        self._gui(self._analysis_finished)

    def _analyse_one_sequence(self, record) -> dict:
        seq   = str(record.seq)
        ssrs  = self.find_ssrs(seq)
        for ssr in ssrs:
            if PRIMER3_AVAILABLE:
                primers = self.design_primers(seq, ssr['start'], ssr['end'])
                ssr['primers'] = primers if primers else {
                    'forward': {'sequence': 'N/A', 'gc_content': 'N/A',
                                'tm': 'N/A', 'size': 'N/A'},
                    'reverse': {'sequence': 'N/A', 'gc_content': 'N/A',
                                'tm': 'N/A', 'size': 'N/A'},
                    'product_size': 'N/A',
                }
            else:
                ssr['primers'] = {
                    'forward': {'sequence': 'N/A', 'gc_content': 'N/A',
                                'tm': 'N/A', 'size': 'N/A'},
                    'reverse': {'sequence': 'N/A', 'gc_content': 'N/A',
                                'tm': 'N/A', 'size': 'N/A'},
                    'product_size': 'N/A',
                }
        return {'id': record.id, 'sequence': seq, 'results': ssrs}

    def design_primers(self, sequence: str, ssr_start: int, ssr_end: int, attempt: int = 1):
        flank_size = 500 if attempt == 1 else 1000
        lf_start   = max(0, ssr_start - flank_size - 1)
        lf_end     = ssr_start - 1
        left_flank = sequence[lf_start:lf_end]
        rf_start   = ssr_end
        rf_end     = min(len(sequence), ssr_end + flank_size)
        right_flank = sequence[rf_start:rf_end]
        ssr_len    = ssr_end - ssr_start + 1

        if len(left_flank) < 50 or len(right_flank) < 50:
            return None

        seq_args = {
            'SEQUENCE_ID':       'SSR',
            'SEQUENCE_TEMPLATE': left_flank + 'N' * ssr_len + right_flank,
            'SEQUENCE_TARGET':   (len(left_flank), ssr_len),
        }

        try:
            params = self.primer_params.copy()
            if attempt == 2:
                params['PRIMER_MIN_GC'] = 20.0
                params['PRIMER_MAX_GC'] = 80.0
                params['PRIMER_MAX_POLY_X'] = 5
                params['PRIMER_MIN_TM'] = 45.0
                params['PRIMER_MAX_TM'] = 75.0

            result = primer3.design_primers(seq_args, params)

            if result.get('PRIMER_PAIR_NUM_RETURNED', 0) > 0:
                return {
                    'forward': {
                        'sequence':   result['PRIMER_LEFT_0_SEQUENCE'],
                        'gc_content': round(result['PRIMER_LEFT_0_GC_PERCENT'], 1),
                        'tm':         round(result['PRIMER_LEFT_0_TM'], 1),
                        'size':       len(result['PRIMER_LEFT_0_SEQUENCE']),
                    },
                    'reverse': {
                        'sequence':   result['PRIMER_RIGHT_0_SEQUENCE'],
                        'gc_content': round(result['PRIMER_RIGHT_0_GC_PERCENT'], 1),
                        'tm':         round(result['PRIMER_RIGHT_0_TM'], 1),
                        'size':       len(result['PRIMER_RIGHT_0_SEQUENCE']),
                    },
                    'product_size': result['PRIMER_PAIR_0_PRODUCT_SIZE'],
                }
            elif attempt == 1:
                return self.design_primers(sequence, ssr_start, ssr_end, attempt=2)
        except Exception as exc:
            logging.warning("Primer3 error at %d-%d (attempt %d): %s",
                           ssr_start, ssr_end, attempt, exc)

        return None

    def _analysis_finished(self):
        self.display_results()
        self.save_results()
        self.status.config(text="Analysis complete!")
        self.analysis_in_progress = False
        self.run_analysis_btn.config(state='normal')

    def find_ssrs(self, sequence: str) -> list:
        all_ssrs = []

        for rtype, info in self.min_repeats.items():
            if not info['enabled']:
                continue

            motif_len = info['motif_length']
            pattern   = self._ssr_patterns.get(rtype)
            if pattern is None:
                continue

            for match in pattern.finditer(sequence):
                start     = match.start() + 1
                end       = match.end()
                raw_motif = match.group(1).upper()
                motif     = self._canonical_motif(raw_motif)

                all_ssrs.append({
                    'type':     rtype,
                    'start':    start,
                    'end':      end,
                    'motif':    motif,
                    'repeats':  (end - start + 1) // motif_len,
                    'length':   end - start + 1,
                    'motif_len': motif_len,
                })

        if not all_ssrs:
            return []

        all_ssrs.sort(key=lambda x: (x['start'], x['motif_len']))

        accepted: list = []
        occupied_end   = 0

        from_start: dict = {}
        for ssr in all_ssrs:
            s = ssr['start']
            if s not in from_start:
                from_start[s] = ssr

        candidates = sorted(from_start.values(), key=lambda x: x['start'])

        occupied_end = 0
        for ssr in candidates:
            if ssr['start'] > occupied_end:
                accepted.append(ssr)
                occupied_end = ssr['end']
            elif ssr['end'] > occupied_end:
                pass

        return accepted

    def _filter_duplicate_primers(self, excel_data: list) -> list:
        seen_pairs = set()
        filtered = []
        for item in excel_data:
            fwd = item.get('Forward Primer', 'N/A')
            rev = item.get('Reverse Primer', 'N/A')
            pair_key = (fwd, rev)
            if pair_key not in seen_pairs and fwd != 'N/A' and rev != 'N/A':
                seen_pairs.add(pair_key)
                filtered.append(item)
        return filtered

    def display_results(self):
        for widget in self.results_frame.winfo_children():
            widget.destroy()

        if not self.analysis_results:
            ttk.Label(self.results_frame, text="No results").pack()
            return

        notebook = ttk.Notebook(self.results_frame)
        notebook.pack(fill=tk.BOTH, expand=True)

        for result in self.analysis_results:
            outer = ttk.Frame(notebook)
            notebook.add(outer, text=result['id'][:20])

            filter_bar = ttk.Frame(outer)
            filter_bar.pack(fill=tk.X, padx=5, pady=(5, 0))
            ttk.Label(filter_bar, text="Filter type:").pack(
                side=tk.LEFT, padx=(0, 4))
            filter_var = tk.StringVar()
            type_choices = ['All', 'mono', 'di', 'tri', 'tetra', 'penta', 'hexa']
            ttk.Combobox(filter_bar, textvariable=filter_var,
                         values=type_choices, state='readonly',
                         width=8).pack(side=tk.LEFT)
            filter_var.set('All')

            ttk.Label(filter_bar, text="  Search:").pack(
                side=tk.LEFT, padx=(10, 4))
            search_var = tk.StringVar()
            ttk.Entry(filter_bar, textvariable=search_var,
                      width=20).pack(side=tk.LEFT)

            ttk.Button(
                filter_bar, text="Export Flanking FASTA",
                command=lambda r=result: self.export_flanking_fasta(r)
            ).pack(side=tk.RIGHT, padx=5)

            columns = [
                'Type', 'Start', 'End', 'Motif', 'Repeats', 'Length',
                'Fwd Primer', 'Fwd GC%', 'Fwd Tm', 'Fwd Size',
                'Rev Primer', 'Rev GC%', 'Rev Tm', 'Rev Size',
                'Product Size', 'Tm Diff',
            ]

            frame = ttk.Frame(outer)
            frame.pack(fill=tk.BOTH, expand=True)

            canvas   = tk.Canvas(frame, borderwidth=0)
            scroll_y = ttk.Scrollbar(frame, orient=tk.VERTICAL,
                                     command=canvas.yview)
            scroll_x = ttk.Scrollbar(frame, orient=tk.HORIZONTAL,
                                     command=canvas.xview)
            canvas.configure(yscrollcommand=scroll_y.set,
                             xscrollcommand=scroll_x.set)

            scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
            scroll_x.pack(side=tk.BOTTOM, fill=tk.X)
            canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

            inner = ttk.Frame(canvas)
            canvas.create_window((0, 0), window=inner, anchor="nw")

            def on_frame_configure(event, c=canvas):
                c.configure(scrollregion=c.bbox("all"))
            inner.bind("<Configure>", on_frame_configure)

            tree = ttk.Treeview(inner, columns=columns,
                                show='headings', height=15)

            for col in columns:
                tree.heading(col, text=col,
                             command=lambda c=col, t=tree: self._sort_tree(t, c))
                tree.column(col, width=100, anchor=tk.CENTER)

            tree.column('Type',       width=50)
            tree.column('Motif',      width=70)
            tree.column('Fwd Primer', width=130)
            tree.column('Rev Primer', width=130)
            tree.pack(fill=tk.BOTH, expand=True)

            all_rows = []
            for ssr in result['results']:
                fwd = ssr.get('primers', {}).get('forward', {})
                rev = ssr.get('primers', {}).get('reverse', {})
                tm_diff = 'N/A'
                if fwd.get('tm', 'N/A') != 'N/A' and \
                        rev.get('tm', 'N/A') != 'N/A':
                    try:
                        tm_diff = f"{abs(float(fwd['tm']) - float(rev['tm'])):.1f}"
                    except Exception:
                        pass
                row = (
                    ssr['type'], ssr['start'], ssr['end'],
                    ssr['motif'], ssr['repeats'], ssr['length'],
                    fwd.get('sequence', 'N/A'),
                    fwd.get('gc_content', 'N/A'),
                    fwd.get('tm', 'N/A'),
                    fwd.get('size', 'N/A'),
                    rev.get('sequence', 'N/A'),
                    rev.get('gc_content', 'N/A'),
                    rev.get('tm', 'N/A'),
                    rev.get('size', 'N/A'),
                    ssr.get('primers', {}).get('product_size', 'N/A'),
                    tm_diff,
                )
                all_rows.append(row)
                tree.insert('', 'end', values=row)

            def _apply_filter(fv=filter_var, sv=search_var,
                              t=tree, rows=all_rows):
                for item in t.get_children():
                    t.delete(item)
                ft = fv.get()
                st = sv.get().lower()
                for row in rows:
                    if ft != 'All' and row[0] != ft:
                        continue
                    if st and not any(
                            st in str(v).lower() for v in row):
                        continue
                    t.insert('', 'end', values=row)

            filter_var.trace_add('write', lambda *_: _apply_filter())
            search_var.trace_add('write', lambda *_: _apply_filter())

            if result['results']:
                self._add_density_tab(notebook, result)

    def _sort_tree(self, tree: ttk.Treeview, col: str):
        data = [(tree.set(item, col), item)
                for item in tree.get_children('')]
        try:
            data.sort(key=lambda x: float(x[0]))
        except ValueError:
            data.sort(key=lambda x: x[0].lower()
                      if isinstance(x[0], str) else x[0])
        for idx, (_, item) in enumerate(data):
            tree.move(item, '', idx)

    def _add_density_tab(self, notebook: ttk.Notebook, result: dict):
        plot_frame = ttk.Frame(notebook)
        notebook.add(plot_frame, text=f"{result['id'][:15]} — Map")

        fig, ax = plt.subplots(figsize=(10, 2.5))
        seq_len  = len(result['sequence'])
        colours  = {'mono': '#e74c3c', 'di': '#3498db', 'tri': '#2ecc71',
                    'tetra': '#f39c12', 'penta': '#9b59b6', 'hexa': '#1abc9c'}

        ax.hlines(0, 0, seq_len, colors='#95a5a6', linewidth=4)
        for ssr in result['results']:
            c = colours.get(ssr['type'], '#34495e')
            ax.vlines(ssr['start'], -0.5, 0.5, colors=c, linewidth=1.5,
                      alpha=0.7)

        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color=colours[t], lw=2, label=t)
            for t in colours if any(
                s['type'] == t for s in result['results'])
        ]
        ax.legend(handles=legend_elements, loc='upper right',
                  fontsize=8, ncol=3)
        ax.set_xlim(0, seq_len)
        ax.set_xlabel("Position (bp)")
        ax.set_yticks([])
        ax.set_title(f"SSR map — {result['id']}", fontsize=10)
        fig.tight_layout()

        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        canvas_widget = FigureCanvasTkAgg(fig, master=plot_frame)
        canvas_widget.draw()
        canvas_widget.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        plt.close(fig)

    def export_flanking_fasta(self, result: dict):
        filepath = filedialog.asksaveasfilename(
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")],
            initialfile=f"{result['id']}_flanking.fasta",
        )
        if not filepath:
            return

        flank = 200
        seq   = result['sequence']
        with open(filepath, 'w') as fh:
            for ssr in result['results']:
                s = max(0, ssr['start'] - 1 - flank)
                e = min(len(seq), ssr['end'] + flank)
                fh.write(
                    f">{result['id']}_{ssr['start']}_{ssr['end']}_"
                    f"{ssr['motif']}x{ssr['repeats']}\n"
                )
                fh.write(seq[s:e] + "\n")

        messagebox.showinfo("Done",
                            f"Flanking FASTA saved to:\n{filepath}")
        logging.info("Flanking FASTA exported: %s", filepath)

    def save_results(self):
        if not self.analysis_results:
            return

        try:
            input_path  = self.file_path.get()
            output_dir  = os.path.dirname(input_path) if input_path else os.getcwd()

            proj        = self.project_name
            if not proj:
                proj = (os.path.splitext(os.path.basename(input_path))[0]
                        if input_path else "output")

            txt_filename   = os.path.join(output_dir, f"{proj}_SSR_results.txt")
            excel_filename = os.path.join(output_dir, f"{proj}_SSR_results.xlsx")

            with open(txt_filename, 'w') as txt_file:
                for result in self.analysis_results:
                    txt_file.write(f"Sequence ID: {result['id']}\n")
                    txt_file.write(f"Sequence Length: {len(result['sequence'])}\n")
                    txt_file.write(f"SSRs Found: {len(result['results'])}\n\n")
                    for ssr in result['results']:
                        primers = ssr.get('primers', {})
                        fwd = primers.get('forward', {})
                        rev = primers.get('reverse', {})
                        tm_diff = 'N/A'
                        if fwd.get('tm', 'N/A') != 'N/A' and \
                                rev.get('tm', 'N/A') != 'N/A':
                            try:
                                tm_diff = (f"{abs(float(fwd['tm']) - float(rev['tm'])):.1f}")
                            except Exception:
                                pass
                        txt_file.write(
                            f"SSR Type: {ssr['type']}\n"
                            f"Position: {ssr['start']}-{ssr['end']}\n"
                            f"Motif: {ssr['motif']} (repeats: {ssr['repeats']}, "
                            f"length: {ssr['length']}bp)\n"
                            f"Forward Primer: {fwd.get('sequence','N/A')}\n"
                            f"  GC%: {fwd.get('gc_content','N/A')}, "
                            f"Tm: {fwd.get('tm','N/A')}°C, "
                            f"Size: {fwd.get('size','N/A')}bp\n"
                            f"Reverse Primer: {rev.get('sequence','N/A')}\n"
                            f"  GC%: {rev.get('gc_content','N/A')}, "
                            f"Tm: {rev.get('tm','N/A')}°C, "
                            f"Size: {rev.get('size','N/A')}bp\n"
                            f"Product Size: {primers.get('product_size','N/A')}bp\n"
                            f"Tm Difference: {tm_diff}°C\n\n"
                        )

            excel_data = []
            marker_counter = 1
            self.marker_name_map = {}

            for result in self.analysis_results:
                for ssr in result['results']:
                    primers = ssr.get('primers', {})
                    fwd = primers.get('forward', {})
                    rev = primers.get('reverse', {})
                    tm_diff_val = 'N/A'
                    if fwd.get('tm', 'N/A') != 'N/A' and \
                            rev.get('tm', 'N/A') != 'N/A':
                        try:
                            tm_diff_val = abs(float(fwd['tm']) - float(rev['tm']))
                        except Exception:
                            pass

                    marker_id = f"Mk{marker_counter}"
                    detail = f"{result['id']}_{ssr['start']}_{ssr['end']}_{ssr['motif']}x{ssr['repeats']}"
                    self.marker_name_map[detail] = marker_id

                    excel_data.append({
                        'Marker':        marker_id,
                        'Details':       detail,
                        'Sequence ID':   result['id'],
                        'SSR Type':      ssr['type'],
                        'Start':         ssr['start'],
                        'End':           ssr['end'],
                        'Motif':         ssr['motif'],
                        'Repeats':       ssr['repeats'],
                        'Length':        ssr['length'],
                        'Forward Primer':fwd.get('sequence', 'N/A'),
                        'Forward GC%':   fwd.get('gc_content', 'N/A'),
                        'Forward Tm':    fwd.get('tm', 'N/A'),
                        'Forward Size':  fwd.get('size', 'N/A'),
                        'Reverse Primer':rev.get('sequence', 'N/A'),
                        'Reverse GC%':   rev.get('gc_content', 'N/A'),
                        'Reverse Tm':    rev.get('tm', 'N/A'),
                        'Reverse Size':  rev.get('size', 'N/A'),
                        'Product Size':  primers.get('product_size', 'N/A'),
                        'Tm Difference': tm_diff_val,
                    })
                    marker_counter += 1

            with pd.ExcelWriter(excel_filename, engine='openpyxl') as writer:
                df_main = pd.DataFrame(excel_data)
                df_main.to_excel(writer,
                                 sheet_name='prospecting and primer design',
                                 index=False)
                self.create_total_summary_sheet(df_main, writer)
                self.create_scaffold_summary_sheet(df_main, writer)

            self.generate_primer_file(excel_data, proj, output_dir)
            logging.info("Results saved: %s, %s", txt_filename, excel_filename)
            self._gui(lambda: self.status.config(
                text=f"Results saved:\n{txt_filename}\n{excel_filename}"))

        except Exception as save_exc:
            logging.error("Failed to save results: %s", save_exc)
            error_msg = str(save_exc)
            self._gui(lambda: messagebox.showerror(
                "Error", f"Failed to save results:\n{error_msg}"))

    def generate_primer_file(self, excel_data: list, base_name: str,
                             output_dir: str):
        try:
            valid_primers = [
                item for item in excel_data
                if item['Forward Primer'] not in ('N/A', '') and
                   item['Reverse Primer'] not in ('N/A', '')
            ]

            valid_primers = self._filter_duplicate_primers(valid_primers)

            primer_data = []
            for item in valid_primers:
                primer_data.append({
                    'Marker':          item['Marker'],
                    'Forward primer':  item['Forward Primer'],
                    'Reverse primer':  item['Reverse Primer'],
                    'Details':         item.get('Details', ''),
                })

            primer_df       = pd.DataFrame(primer_data)
            primer_filename = os.path.join(output_dir, f"{base_name}_Primers.xlsx")
            primer_df.to_excel(primer_filename, index=False)

            wb = openpyxl.load_workbook(primer_filename)
            ws = wb.active
            for col in ws.columns:
                max_len = max(
                    (len(str(cell.value or '')) for cell in col), default=0)
                ws.column_dimensions[
                    col[0].column_letter].width = min(max_len + 2, 50)
            for cell in ws[1]:
                cell.font = Font(bold=True)
            wb.save(primer_filename)

        except Exception as gen_exc:
            logging.warning("Could not create primer file: %s", gen_exc)
            error_msg = str(gen_exc)
            self._gui(lambda: messagebox.showwarning("Warning",
                                   f"Could not create primer file: {error_msg}"))

    def create_total_summary_sheet(self, df: pd.DataFrame, writer):
        ssr_types = ['mono', 'di', 'tri', 'tetra', 'penta', 'hexa']
        total_counts   = {t: len(df[df['SSR Type'] == t]) for t in ssr_types}
        primer_counts  = {
            t: len(df[(df['SSR Type'] == t) &
                      (df['Forward Primer'] != 'N/A') &
                      (df['Reverse Primer'] != 'N/A')])
            for t in ssr_types
        }
        summary = pd.DataFrame({
            'SSR Type':               ['Mono', 'Di', 'Tri', 'Tetra', 'Penta', 'Hexa'],
            'Total SSRs Identified':  [total_counts[t]  for t in ssr_types],
            'SSRs with Primers':      [primer_counts[t] for t in ssr_types],
        })
        summary.to_excel(writer, sheet_name='Total summary', index=False)
        ws = writer.sheets['Total summary']
        for ref_col, title, anchor in [
            (2, 'Total SSRs Identified', 'E2'),
            (3, 'SSRs with Primers',     'E20'),
        ]:
            chart = BarChart()
            chart.type  = "col"
            chart.style = 10
            chart.title = title
            chart.y_axis.title = "Number of SSRs"
            chart.x_axis.title = "SSR Type"
            chart.legend = None
            data = Reference(ws, min_col=ref_col,
                             min_row=1, max_row=7)
            cats = Reference(ws, min_col=1, min_row=2, max_row=7)
            chart.add_data(data, titles_from_data=True)
            chart.set_categories(cats)
            chart.width  = 15
            chart.height = 10
            ws.add_chart(chart, anchor)

    def create_scaffold_summary_sheet(self, df: pd.DataFrame, writer):
        if 'Sequence ID' not in df.columns:
            return
        rows = []
        for scaffold in df['Sequence ID'].unique():
            sub = df[df['Sequence ID'] == scaffold]
            for t in ['mono', 'di', 'tri', 'tetra', 'penta', 'hexa']:
                total  = len(sub[sub['SSR Type'] == t])
                with_p = len(sub[
                    (sub['SSR Type'] == t) &
                    (sub['Forward Primer'] != 'N/A') &
                    (sub['Reverse Primer'] != 'N/A')
                ])
                rows.append({
                    'Scaffold':                 scaffold,
                    'SSR Type':                 t,
                    'Total SSRs Identified':    total,
                    'SSRs with Primers':        with_p,
                })
            rows.append({'Scaffold': '', 'SSR Type': '',
                         'Total SSRs Identified': '',
                         'SSRs with Primers': ''})

        pd.DataFrame(rows).to_excel(
            writer, sheet_name='Summary by scaffold', index=False)
        ws  = writer.sheets['Summary by scaffold']
        row = 2
        for scaffold in df['Sequence ID'].unique():
            count = sum(1 for r in rows if r['Scaffold'] == scaffold)
            for c in range(1, 5):
                ws.cell(row=row, column=c).font = Font(bold=True)
            row += count + 1

    # ==================== MÉTODOS DE BLAST ====================

    def browse_blast_input(self):
        fp = filedialog.askopenfilename(
            filetypes=[('Excel files', '*.xlsx *.xls'), ('All files', '*.*')])
        if fp:
            self.blast_input_var.set(fp)
            self.blast_input_file = fp
            self.current_results_df = pd.read_excel(fp)

    def _get_original_name_from_marker(self, marker_id: str, df: pd.DataFrame) -> str:
        if 'Details' in df.columns:
            marker_rows = df[df['Marker'] == marker_id]
            if not marker_rows.empty:
                details = marker_rows.iloc[0]['Details']
                if pd.notna(details):
                    return str(details)
        if marker_id in self.marker_name_map.values():
            for detail, mk in self.marker_name_map.items():
                if mk == marker_id:
                    return detail
        return marker_id

    def _get_seq_id_from_original_name(self, original_name: str) -> str:
        if not original_name or original_name == 'N/A':
            return None

        if original_name in self.seq_id_map:
            return self.seq_id_map[original_name]

        if self.id_mapping and original_name in self.id_mapping:
            return self.id_mapping[original_name]

        pattern = r'^(.+?)(?:_\d+_+\d+_+[ACGT]+x\d+)$'
        match = re.match(pattern, original_name)
        if match:
            candidate = match.group(1)
            if candidate in self.reference_sequences:
                return candidate
            if candidate in self.seq_id_map:
                return self.seq_id_map[candidate]

        parts = original_name.split('_')
        for i, part in enumerate(parts):
            if part.isdigit() and i > 0:
                seq_id = '_'.join(parts[:i])
                if seq_id in self.reference_sequences:
                    return seq_id
                for ref_id in self.reference_sequences:
                    if seq_id in ref_id or ref_id in seq_id:
                        return ref_id
                break

        for ref_id in self.reference_sequences:
            if original_name.startswith(ref_id) or ref_id in original_name:
                return ref_id

        return parts[0] if parts else None

    def check_gene_link(self, seq_id: str, marker_start: int,
                        marker_end: int, gene_flank: int,
                        min_overlap: int) -> str:
        if seq_id not in self.blast_genbank_features:
            return "No annotation available"
        linked = []
        for feat in self.blast_genbank_features[seq_id]:
            exp_s = max(1, feat['start'] - gene_flank)
            exp_e = feat['end'] + gene_flank
            ov_s  = max(marker_start, exp_s)
            ov_e  = min(marker_end,   exp_e)
            if ov_e > ov_s:
                ov = ov_e - ov_s + 1
                if ov >= min_overlap:
                    linked.append(
                        f"{feat['gene']} ({feat['type']}, "
                        f"strand {feat['strand']}, overlap {ov}bp)"
                    )
        return " | ".join(linked) if linked else "Intergenic region"

    def _load_blast_df(self):
        if self.current_results_df is not None:
            df = self.current_results_df
            required = {'Marker', 'Start', 'End'}
            if not required.issubset(df.columns):
                required_alt = {'Sequence ID', 'Start', 'End'}
                if required_alt.issubset(df.columns):
                    df = df.rename(columns={'Sequence ID': 'Marker'})
                else:
                    messagebox.showerror(
                        "Error",
                        "Excel file must contain 'Marker', 'Start', and 'End' columns")
                    return None
            return df

        if not self.blast_input_file:
            messagebox.showwarning(
                "Warning",
                "Please select an Excel file with SSRs and primers")
            return None
        try:
            df = pd.read_excel(self.blast_input_file)
            self.current_results_df = df
            required = {'Marker', 'Start', 'End'}
            if not required.issubset(df.columns):
                required_alt = {'Sequence ID', 'Start', 'End'}
                if required_alt.issubset(df.columns):
                    df = df.rename(columns={'Sequence ID': 'Marker'})
                else:
                    messagebox.showerror(
                        "Error",
                        "Excel file must contain 'Marker', 'Start', and 'End' columns")
                    return None
            return df
        except Exception as load_exc:
            error_msg = str(load_exc)
            messagebox.showerror("Error",
                                 f"Failed to read Excel file:\n{error_msg}")
            return None

    def run_blast_analysis(self):
        if not self.family_var.get():
            messagebox.showwarning("Warning",
                                   "Please select a family for BLAST")
            return
        if not self.blast_fasta_file:
            messagebox.showwarning("Warning",
                                   "Please load the reference FASTA file first")
            return
        if not self.blast_input_file and self.current_results_df is None:
            messagebox.showwarning("Warning",
                                   "Please load the SSR Results file first")
            return

        try:
            flank = int(self.flank_size_var.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid flanking region size")
            return
        df = self._load_blast_df()
        if df is None:
            return

        self.blast_in_progress       = True
        self.blast_progress['value'] = 0
        self.blast_progress['maximum'] = len(df)
        self.blast_status.config(text="Starting NCBI BLAST…")

        database = self.blast_db_var.get()
        if database == "mRNA_refSeq":
            db = "refseq_rna"
        else:
            db = "nt"

        self.run_ncbi_blast_btn.config(state='disabled')
        self.run_local_blast_btn.config(state='disabled')

        threading.Thread(
            target=self.perform_blast_analysis,
            args=(flank, df, True, db), daemon=True
        ).start()

    def run_local_blast_analysis(self):
        if not self.blast_genbank_dict and not self.blast_genbank_features:
            if self.blast_genbank_file:
                self.browse_genbank()
            elif self.blast_gff_file and self.gff_fasta_var.get():
                self.load_gff_annotations(self.blast_gff_file, self.gff_fasta_var.get())
            elif self.blast_gtf_file and self.gtf_fasta_var.get():
                self.load_gtf_annotations(self.blast_gtf_file, self.gtf_fasta_var.get())
            else:
                messagebox.showwarning("Warning",
                                   "Please load a GenBank, GFF3, or GTF file first")
                return
        if not self.blast_fasta_file:
            messagebox.showwarning("Warning",
                                   "Please load the reference FASTA file first")
            return
        if not self.blast_input_file and self.current_results_df is None:
            messagebox.showwarning("Warning",
                                   "Please load the SSR Results file first")
            return
        try:
            flank      = int(self.flank_size_var.get())
            gene_flank = int(self.gene_flank_var.get())
            min_ov     = int(self.min_overlap_var.get())
            max_mm     = float(self.max_mismatch_var.get()) / 100.0
        except ValueError as exc:
            messagebox.showerror("Error", f"Invalid parameter: {exc}")
            return
        df = self._load_blast_df()
        if df is None:
            return

        if not self._kmer_indices_built:
            self._build_all_kmer_indices()

        self.blast_in_progress        = True
        self.blast_progress['value']  = 0
        self.blast_progress['maximum']= len(df)
        self.blast_status.config(text="Starting Local BLAST…")

        self.run_ncbi_blast_btn.config(state='disabled')
        self.run_local_blast_btn.config(state='disabled')

        threading.Thread(
            target=self.perform_local_blast_analysis,
            args=(flank, df, max_mm, gene_flank, min_ov),
            daemon=True
        ).start()

    def perform_blast_analysis(self, flank_size: int, df: pd.DataFrame,
                               use_ncbi: bool, database: str):
        results   = []
        processed = 0

        for _, row in df.iterrows():
            if not self.blast_in_progress:
                break

            marker_id = row['Marker']
            start  = int(row['Start'])
            end    = int(row['End'])

            try:
                original_name = self._get_original_name_from_marker(marker_id, df)
                seq_id = self._get_seq_id_from_original_name(original_name)

                if seq_id is None or seq_id not in self.reference_sequences:
                    results.append({
                        'marker': marker_id,
                        'accession': f'Sequence not found: {seq_id}', 'description': '',
                        'e_value': '', 'identity': '', 'alignment': '',
                        'gene_info': 'N/A', 'original_name': original_name,
                    })
                    processed += 1
                    _p = processed
                    self._gui(lambda p=_p: self.blast_progress.__setitem__(
                        'value', p))
                    continue

                full_seq = self.reference_sequences[seq_id]
                r_start  = max(0, start - flank_size - 1)
                r_end    = min(len(full_seq), end + flank_size)
                region   = full_seq[r_start:r_end]

                _m = marker_id
                self._gui(lambda m=_m: self.blast_status.config(
                    text=f"NCBI BLAST: {m}…"))
                fam = self.family_var.get()
                if ": " in fam:
                    fam = fam.split(": ")[1]

                try:
                    handle = NCBIWWW.qblast(
                        "blastn", database, region,
                        entrez_query=f"{fam}[Organism]",
                        hitlist_size=5, expect=10,
                        word_size=11, nucl_penalty=-3, nucl_reward=2,
                    )
                    best = None
                    for blast_record in NCBIXML.parse(handle):
                        for aln in blast_record.alignments:
                            for hsp in aln.hsps:
                                if best is None or hsp.expect < best['e_value']:
                                    best = {
                                        'accession': aln.accession,
                                        'title':     aln.title,
                                        'e_value':   hsp.expect,
                                        'identity':  hsp.identities / hsp.align_length * 100,
                                        'alignment': str(hsp),
                                    }
                except socket.timeout:
                    logging.warning("NCBI BLAST timeout for %s", marker_id)
                    results.append({
                        'marker': marker_id, 'accession': 'Timeout',
                        'description': 'NCBI BLAST timed out', 'e_value': '',
                        'identity': '', 'alignment': '', 'gene_info': 'N/A',
                        'original_name': original_name,
                    })
                    processed += 1
                    _p = processed
                    self._gui(lambda p=_p: self.blast_progress.__setitem__('value', p))
                    continue
                except Exception as blast_exc:
                    raise blast_exc

                if best:
                    results.append({
                        'marker':      marker_id,
                        'accession':   best['accession'],
                        'description': best['title'],
                        'e_value':     best['e_value'],
                        'identity':    best['identity'],
                        'alignment':   best['alignment'],
                        'gene_info':   'N/A (NCBI BLAST)',
                        'original_name': original_name,
                    })
                else:
                    results.append({
                        'marker': marker_id, 'accession': 'No hits',
                        'description': '', 'e_value': '',
                        'identity': '', 'alignment': '', 'gene_info': 'N/A',
                        'original_name': original_name,
                    })
            except Exception as blast_exc:
                logging.error("NCBI BLAST error for %s: %s", marker_id, blast_exc)
                error_msg = str(blast_exc)
                results.append({
                    'marker': marker_id, 'accession': 'BLAST failed',
                    'description': error_msg, 'e_value': '',
                    'identity': '', 'alignment': '', 'gene_info': 'N/A',
                    'original_name': original_name if 'original_name' in locals() else marker_id,
                })

            processed += 1
            _p = processed
            self._gui(lambda p=_p: self.blast_progress.__setitem__('value', p))
            time.sleep(0.5)

        self._gui(lambda: self.run_ncbi_blast_btn.config(state='normal'))
        self._gui(lambda: self.run_local_blast_btn.config(state='normal'))

        if self.blast_in_progress:
            self.save_blast_results(results, "NCBI")
            self._gui(lambda: self.blast_status.config(
                text="NCBI BLAST analysis complete!"))
            self.blast_in_progress = False

    def perform_local_blast_analysis(self, flank_size: int, df: pd.DataFrame,
                                     max_mismatch: float, gene_flank: int,
                                     min_overlap: int):
        results   = []
        processed = 0
        total     = len(df)

        for _, row in df.iterrows():
            if not self.blast_in_progress:
                break

            marker_id = row['Marker']
            start  = int(row['Start'])
            end    = int(row['End'])

            try:
                original_name = self._get_original_name_from_marker(marker_id, df)
                seq_id = self._get_seq_id_from_original_name(original_name)

                if seq_id is None or seq_id not in self.reference_sequences:
                    results.append({
                        'marker': marker_id,
                        'accession': f'Sequence not found: {seq_id}',
                        'description': 'Source sequence not found in FASTA',
                        'e_value': '', 'identity': '', 'alignment': '',
                        'gene_info': 'N/A', 'status': 'ERROR',
                        'original_name': original_name,
                    })
                    processed += 1
                    _p = processed
                    self._gui(lambda p=_p: self.blast_progress.__setitem__(
                        'value', p))
                    continue

                full_seq = self.reference_sequences[seq_id]
                r_start  = max(0, start - flank_size - 1)
                r_end    = min(len(full_seq), end + flank_size)
                region   = full_seq[r_start:r_end]
                region_rev = self.reverse_complement(region)

                _m, _t = marker_id, total
                _p2 = processed + 1
                self._gui(lambda m=_m, p=_p2, t=_t: self.blast_status.config(
                    text=f"Local BLAST: {m} ({p}/{t})…"))

                all_matches = []

                for gb_id, gb_seq in self.blast_genbank_dict.items():
                    hits_fwd = self._find_with_mismatch(
                        region, gb_id, gb_seq, max_mismatch)
                    hits_rev = self._find_with_mismatch(
                        region_rev, gb_id, gb_seq, max_mismatch)

                    for pos, identity in hits_fwd:
                        gene_info = self.check_gene_link(
                            gb_id, pos, pos + len(region) - 1,
                            gene_flank, min_overlap)
                        all_matches.append({
                            'accession': gb_id,
                            'position':  pos,
                            'identity':  identity,
                            'alignment': f"Match at position {pos} with "
                                         f"{identity:.1f}% identity (forward)",
                            'gene_info': gene_info,
                        })

                    for pos, identity in hits_rev:
                        gene_info = self.check_gene_link(
                            gb_id, pos, pos + len(region_rev) - 1,
                            gene_flank, min_overlap)
                        all_matches.append({
                            'accession': gb_id,
                            'position':  pos,
                            'identity':  identity,
                            'alignment': f"Match at position {pos} with "
                                         f"{identity:.1f}% identity (reverse complement)",
                            'gene_info': gene_info,
                        })

                all_matches.sort(key=lambda x: -x['identity'])

                if all_matches:
                    best = all_matches[0]
                    results.append({
                        'marker':        marker_id,
                        'accession':     best['accession'],
                        'description':   f"Position: {best['position']}",
                        'e_value':       'N/A',
                        'identity':      best['identity'],
                        'alignment':     best['alignment'],
                        'gene_info':     best['gene_info'],
                        'total_matches': len(all_matches),
                        'status':        'FOUND',
                        'original_name': original_name,
                    })
                else:
                    results.append({
                        'marker':    marker_id,
                        'accession': 'No match found',
                        'description': f'No match ≥{(1-max_mismatch)*100:.0f}%',
                        'e_value': '', 'identity': '', 'alignment': '',
                        'gene_info': 'N/A', 'status': 'NOT_FOUND',
                        'original_name': original_name,
                    })
            except Exception as local_exc:
                logging.error("Local BLAST error for %s: %s", marker_id, local_exc)
                error_msg = str(local_exc)
                results.append({
                    'marker': marker_id,
                    'accession': 'BLAST failed',
                    'description': error_msg, 'e_value': '',
                    'identity': '', 'alignment': '', 'gene_info': 'N/A',
                    'status': 'ERROR',
                    'original_name': original_name if 'original_name' in locals() else marker_id,
                })

            processed += 1
            _p = processed
            self._gui(lambda p=_p: self.blast_progress.__setitem__('value', p))

        self._gui(lambda: self.run_ncbi_blast_btn.config(state='normal'))
        self._gui(lambda: self.run_local_blast_btn.config(state='normal'))

        if self.blast_in_progress:
            self.save_blast_results(results, "Local")
            self._gui(lambda: self.blast_status.config(
                text="Local BLAST complete!"))
            self.blast_in_progress = False

    def cancel_blast(self):
        self.blast_in_progress = False
        self.blast_status.config(text="BLAST analysis cancelled")
        self.run_ncbi_blast_btn.config(state='normal')
        self.run_local_blast_btn.config(state='normal')
        logging.info("BLAST cancelled by user")

    def save_blast_results(self, blast_results: list, blast_type: str):
        if not blast_results:
            return
        try:
            base   = os.path.splitext(os.path.basename(
                self.blast_input_file or 'output'))[0]
            outdir = os.path.dirname(
                self.blast_input_file) if self.blast_input_file else os.getcwd()
            ts     = time.strftime("%Y%m%d_%H%M%S")
            fname  = os.path.join(outdir,
                                  f"{base}_{blast_type}_BLAST_{ts}.xlsx")

            df_results = pd.DataFrame(blast_results)

            if 'original_name' in df_results.columns:
                marker_col_idx = df_results.columns.get_loc('marker') if 'marker' in df_results.columns else 0
                original_name_col = df_results.pop('original_name')
                insert_pos = min(marker_col_idx + 1, len(df_results.columns))
                df_results.insert(insert_pos, 'Original_Marker_Name', original_name_col)

            df_results.to_excel(fname, index=False)
            
            # ========== REMOVIDO: display_blast_results ==========
            # Os resultados são salvos em Excel, não exibidos na interface
            # self._gui(lambda: self.display_blast_results(blast_results))
            
            self._gui(lambda: self.blast_status.config(
                text=f"BLAST results saved to:\n{fname}"))
            self._gui(lambda: messagebox.showinfo(
                "Success", f"BLAST results saved to:\n{fname}"))
            logging.info("BLAST results saved: %s", fname)
        except Exception as save_exc:
            logging.error("Failed to save BLAST results: %s", save_exc)
            error_msg = str(save_exc)
            self._gui(lambda: messagebox.showerror(
                "Error", f"Failed to save BLAST results:\n{error_msg}"))

    # ==================== MÉTODOS DE VALIDAÇÃO ====================

    def toggle_relationship_options(self):
        if self.compute_relationship_var.get():
            self.relationship_frame.pack(fill=tk.X, pady=5)
        else:
            self.relationship_frame.pack_forget()

    def validate_primers(self):
        if not all([self.fasta_file.get(), self.excel_file.get()]):
            messagebox.showerror("Error", "Please select both input files")
            return

        try:
            self.validation_status.config(text="Loading genomes…")
            self.root.update()

            genomes = {}
            for record in SeqIO.parse(self.fasta_file.get(), "fasta"):
                genomes[record.description.split()[0]] = str(record.seq).upper()
            if not genomes:
                messagebox.showerror("Error",
                                     "No valid genomes found in the FASTA file")
                return

            wb    = openpyxl.load_workbook(self.excel_file.get())
            ws    = wb.active
            primers = []
            for row in ws.iter_rows(min_row=2, values_only=True):
                if row[0] and row[1] and row[2]:
                    primers.append({
                        'locus':   str(row[0]),
                        'f_primer':str(row[1]).upper(),
                        'r_primer':str(row[2]).upper(),
                    })
            if not primers:
                messagebox.showerror("Error",
                                     "No valid primers found in the Excel file")
                return

            output_path = self.output_file.get()
            output_dir  = os.path.dirname(output_path)
            os.makedirs(output_dir, exist_ok=True)

            gel_dir = os.path.join(output_dir, "gel_images")
            if self.generate_gel_var.get():
                os.makedirs(gel_dir, exist_ok=True)

            out_wb = openpyxl.Workbook()
            out_wb.remove(out_wb.active)

            all_amplicons        = defaultdict(list)
            all_species_gel_data = {}
            no_amp_genomes       = []

            if self.compute_relationship_var.get():
                presence_matrix = {s: [] for s in genomes}
                all_alleles     = []

            total   = len(genomes)
            current = 0

            for species, genome_seq in genomes.items():
                current += 1
                self.validation_progress['value'] = (current / total) * 50
                self.validation_status.config(
                    text=f"Analysing {species} ({current}/{total})…")
                self.root.update()

                sws     = out_wb.create_sheet(title=species[:31])
                headers = ["Locus", "Forward Primer",
                           "Reverse Primer", "Amplicon Size (bp)"]
                for ci, h in enumerate(headers, 1):
                    cell = sws[f"{get_column_letter(ci)}1"]
                    cell.value = h
                    cell.font  = Font(bold=True)

                gel_data = {p['locus']: [] for p in primers}
                all_species_gel_data[species] = gel_data

                if self.compute_relationship_var.get():
                    species_alleles = []

                row_num  = 2
                has_amps = False
                for primer in primers:
                    fp_seq  = primer['f_primer']
                    rp_rc   = self.reverse_complement(primer['r_primer'])

                    f_pos = [i for i in range(len(genome_seq))
                             if genome_seq.startswith(fp_seq, i)]
                    r_pos = [i for i in range(len(genome_seq))
                             if genome_seq.startswith(rp_rc, i)]

                    amp_sizes = sorted([
                        rp - fp + len(primer['r_primer'])
                        for fp in f_pos for rp in r_pos if rp > fp
                    ])

                    sws[f"A{row_num}"] = primer['locus']
                    sws[f"B{row_num}"] = primer['f_primer']
                    sws[f"C{row_num}"] = primer['r_primer']

                    if amp_sizes:
                        sws[f"D{row_num}"] = ", ".join(map(str, amp_sizes))
                        gel_data[primer['locus']] = amp_sizes
                        has_amps = True
                        all_amplicons[primer['locus']].extend(amp_sizes)
                        if self.compute_relationship_var.get():
                            aid = f"{primer['locus']}_{min(amp_sizes)}"
                            species_alleles.append(aid)
                            if aid not in all_alleles:
                                all_alleles.append(aid)
                    else:
                        sws[f"D{row_num}"] = "Not found"
                        if self.compute_relationship_var.get():
                            aid = f"{primer['locus']}_0"
                            species_alleles.append(aid)
                            if aid not in all_alleles:
                                all_alleles.append(aid)

                    row_num += 1

                if not has_amps:
                    no_amp_genomes.append(species)
                    out_wb.remove(sws)
                    del all_species_gel_data[species]
                else:
                    if self.compute_relationship_var.get():
                        presence_matrix[species] = species_alleles
                    for col in sws.columns:
                        ml  = max((len(str(c.value or '')) for c in col),
                                  default=0)
                        sws.column_dimensions[
                            col[0].column_letter].width = (ml + 2) * 1.2

            if (self.generate_gel_var.get() and
                    self.gel_by_scaffold_var.get() and
                    all_species_gel_data):
                self.validation_status.config(text="Generating scaffold gels…")
                self.root.update()
                self.generate_scaffold_gels(
                    all_species_gel_data, primers, gel_dir)

            if (self.generate_gel_var.get() and
                    self.gel_by_marker_var.get() and
                    all_species_gel_data):
                self.validation_status.config(text="Generating marker gels…")
                self.root.update()
                self.generate_marker_gels(
                    all_species_gel_data, primers, gel_dir)

            sum_ws = out_wb.create_sheet(title="Amplicon_Summary", index=0)
            sum_ws.append(["Marker", "Forward Primer",
                           "Reverse Primer", "Amplicon Sizes (bp)"])
            for primer in primers:
                locus = primer['locus']
                if locus in all_amplicons:
                    sizes = sorted(set(all_amplicons[locus]))
                    sum_ws.append([locus, primer['f_primer'],
                                   primer['r_primer'],
                                   ", ".join(map(str, sizes))])
            for cell in sum_ws[1]:
                cell.font = Font(bold=True)
            for col in sum_ws.columns:
                ml = max((len(str(cell.value or '')) for cell in col), default=0)
                sum_ws.column_dimensions[
                    col[0].column_letter].width = (ml + 2) * 1.2

            if no_amp_genomes:
                na_ws = out_wb.create_sheet(title="No_Amplicons")
                na_ws.append(["Genomes with no amplicons:"])
                for g in no_amp_genomes:
                    na_ws.append([g])
                na_ws['A1'].font = Font(bold=True)
                na_ws.column_dimensions['A'].width = 50

            if self.compute_relationship_var.get():
                self.validation_status.config(
                    text="Computing genetic relationships…")
                self.root.update()

                species_names = list(genomes.keys())
                def _aid_locus(a): return a.rsplit('_', 1)[0]
                def _aid_size(a):  return int(a.rsplit('_', 1)[1])

                loci         = sorted({_aid_locus(a) for a in all_alleles})
                all_possible = []
                for locus in loci:
                    sizes = {
                        _aid_size(a)
                        for s in species_names
                        for a in presence_matrix.get(s, [])
                        if _aid_locus(a) == locus and _aid_size(a) > 0
                    }
                    for sz in sorted(sizes):
                        all_possible.append(f"{locus}_{sz}")
                    all_possible.append(f"{locus}_0")

                binary_data = [
                    [1 if a in presence_matrix.get(s, []) else 0
                     for a in all_possible]
                    for s in species_names
                ]
                binary_df = pd.DataFrame(binary_data,
                                         index=species_names,
                                         columns=all_possible)

                pa_ws = out_wb.create_sheet(title="Presence_Absence")
                pa_ws.cell(1, 1, "Species")
                for ci, a in enumerate(all_possible, 2):
                    pa_ws.cell(1, ci, a)
                for ri, s in enumerate(species_names, 2):
                    pa_ws.cell(ri, 1, s)
                    for ci, a in enumerate(all_possible, 2):
                        pa_ws.cell(ri, ci,
                                   binary_df.loc[s, a])

                jac_sim = 1 - pairwise_distances(
                    binary_df.values, metric='jaccard')
                jac_df  = pd.DataFrame(jac_sim,
                                       index=species_names,
                                       columns=species_names)

                sm_ws = out_wb.create_sheet(title="Similarity_Matrix")
                sm_ws.cell(1, 1, "Species")
                for ci, s in enumerate(species_names, 2):
                    sm_ws.cell(1, ci, s)
                for ri, s1 in enumerate(species_names, 2):
                    sm_ws.cell(ri, 1, s1)
                    for ci, s2 in enumerate(species_names, 2):
                        sm_ws.cell(ri, ci, jac_df.loc[s1, s2])

                self.validation_progress['value'] = 75
                self.validation_status.config(
                    text="Generating dendrogram…")
                self.root.update()

                dend_dir = os.path.join(output_dir, "dendrograms")
                os.makedirs(dend_dir, exist_ok=True)
                jac_dist = 1 - jac_df

                bs_vals = None
                if self.bootstrap_var.get():
                    self.validation_status.config(
                        text="Bootstrapping…")
                    self.root.update()
                    bs_vals = self.perform_bootstrapping(
                        binary_df, reps=self.bootstrap_reps.get())

                self.generate_dendrogram(
                    jac_dist, species_names, dend_dir, bs_vals)
                self.validation_progress['value'] = 90

            try:
                out_wb.save(output_path)
                self.validation_progress['value'] = 100
                self.validation_status.config(
                    text=f"Validation complete! Saved to {output_path}")
                messagebox.showinfo(
                    "Success", "Primer validation completed successfully!")
                logging.info("Validation complete: %s", output_path)
            except PermissionError:
                messagebox.showerror(
                    "Error",
                    "Permission denied. Close the Excel file if it is open.")
                self.validation_status.config(text="Error: could not save")

        except Exception as val_exc:
            logging.error("Validation failed: %s", val_exc)
            error_msg = str(val_exc)
            messagebox.showerror("Error", f"Validation failed: {error_msg}")
            self.validation_status.config(text="Error during validation")

    def generate_scaffold_gels(self, all_species_data: dict,
                               all_primers: list, output_dir: str):
        LANES = 10
        for species, gel_data in all_species_data.items():
            chunks = [all_primers[i:i + LANES]
                      for i in range(0, len(all_primers), LANES)]
            for gel_num, chunk in enumerate(chunks, 1):
                while len(chunk) < LANES:
                    chunk.append({'locus': 'empty'})
                self._draw_gel(
                    chunk, gel_data, 'red',
                    f"Scaffold Gel — {species}" + (
                        f" (Part {gel_num}/{len(chunks)})"
                        if len(chunks) > 1 else ""),
                    os.path.join(
                        output_dir,
                        f"{''.join(c if c.isalnum() else '_' for c in species)}"
                        f"_scaffold_gel_{gel_num}.jpg"),
                )

    def generate_marker_gels(self, all_species_data: dict,
                             all_primers: list, output_dir: str):
        LANES       = 10
        all_species = list(all_species_data.keys())
        for primer in all_primers:
            locus      = primer['locus']
            marker_dat = {s: all_species_data[s].get(locus, [])
                          for s in all_species}
            chunks     = [all_species[i:i + LANES]
                          for i in range(0, len(all_species), LANES)]
            for gel_num, chunk in enumerate(chunks, 1):
                gel_data_chunk = {s: marker_dat.get(s, []) for s in chunk}
                while len(chunk) < LANES:
                    chunk.append('empty')
                self._draw_gel(
                    [{'locus': s} for s in chunk],
                    gel_data_chunk, 'blue',
                    f"Marker Gel — {locus}" + (
                        f" (Part {gel_num}/{len(chunks)})"
                        if len(chunks) > 1 else ""),
                    os.path.join(
                        output_dir,
                        f"marker_{''.join(c if c.isalnum() else '_' for c in locus)}"
                        f"_gel_{gel_num}.jpg"),
                )

    def _draw_gel(self, lanes_info: list, gel_data: dict,
                  band_colour: str, title: str, output_path: str):
        ladder    = [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100]
        num_lanes = len(lanes_info) + 1
        fig, ax   = plt.subplots(figsize=(12, 8))
        plt.subplots_adjust(top=0.85)

        ax.add_patch(Rectangle((0, 0), num_lanes, 1100,
                               facecolor='#D0E0F0',
                               edgecolor='black', alpha=0.7))

        for i in range(num_lanes):
            ax.add_patch(Rectangle((i + 0.1, 1050), 0.8, 30,
                                   facecolor='white',
                                   edgecolor='black', lw=0.5))

        for sz in ladder:
            ax.add_patch(Rectangle((0.2, sz), 0.6, 8,
                                   facecolor='black', edgecolor='black',
                                   lw=0.5))
            ax.text(-0.05, sz, str(sz), ha='right', va='center',
                    fontsize=8, color='black')

        has_large = False
        for i, lane in enumerate(lanes_info):
            lx    = i + 1
            locus = lane.get('locus', 'empty')
            if locus == 'empty':
                ax.text(lx + 0.5, 1100, "Empty", ha='center',
                        va='bottom', fontsize=8, color='gray')
                continue

            sizes = gel_data.get(locus, [])
            ax.text(lx + 0.5, 1100, locus, ha='center', va='bottom',
                    fontsize=8, color='black',
                    bbox=dict(facecolor='white', edgecolor='none',
                              pad=1, alpha=0.7))

            clipped = []
            for sz in sizes:
                if sz > 1000:
                    has_large = True
                    clipped.append(1000)
                else:
                    clipped.append(sz)

            for sz in sorted(set(clipped)):
                ax.add_patch(Rectangle((lx + 0.2, sz), 0.6, 8,
                                       facecolor=band_colour,
                                       edgecolor='black', lw=0.5))

        if has_large:
            ax.text(num_lanes - 0.2, 1020,
                    "* = amplicons >1000 bp",
                    ha='right', va='bottom', fontsize=8,
                    color=band_colour,
                    bbox=dict(facecolor='white',
                              edgecolor=band_colour, pad=2))

        ax.set_xlim(0, num_lanes)
        ax.set_ylim(0, 1150)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        plt.title(title, fontsize=12, fontweight='bold', pad=15)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

    def perform_bootstrapping(self, binary_df: pd.DataFrame,
                              reps: int = 100) -> np.ndarray:
        X        = binary_df.values
        n        = len(X)
        orig_d   = pairwise_distances(X, metric='jaccard')
        orig_Z   = linkage(orig_d, method='average')
        counts   = np.zeros(n - 1)

        for _ in range(reps):
            idx   = np.random.choice(X.shape[1],
                                     size=X.shape[1], replace=True)
            Xr    = X[:, idx]
            rd    = pairwise_distances(Xr, metric='jaccard')
            rZ    = linkage(rd, method='average')

            orig_a = cut_tree(orig_Z, n_clusters=np.arange(2, n + 1))
            res_a  = cut_tree(rZ,    n_clusters=np.arange(2, n + 1))

            for j in range(n - 1):
                if j < orig_a.shape[1] and j < res_a.shape[1]:
                    if np.array_equal(orig_a[:, j], res_a[:, j]):
                        counts[j] += 1

        return (counts / reps) * 100

    def generate_dendrogram(self, distance_matrix, species_names: list,
                            output_dir: str,
                            bootstrap_values=None):
        try:
            os.makedirs(output_dir, exist_ok=True)
            labels = [n.replace('_', ' ') for n in species_names]
            Z      = linkage(distance_matrix.values, method='average',
                             optimal_ordering=True)

            fig, ax = plt.subplots(figsize=(12, 8))
            for spine in ax.spines.values():
                spine.set_visible(False)
            plt.title('UPGMA Dendrogram', fontsize=14, pad=20)
            plt.xlabel('Species', fontsize=12)
            plt.ylabel('Genetic Distance', fontsize=12)

            ddata = dendrogram(Z, labels=labels,
                               leaf_rotation=90,
                               leaf_font_size=10,
                               orientation='top',
                               above_threshold_color='grey', ax=ax)

            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center')
            plt.subplots_adjust(bottom=0.2)

            if bootstrap_values is not None:
                icoord = np.array(ddata['icoord'])
                dcoord = np.array(ddata['dcoord'])
                for i, (ic, dc) in enumerate(zip(icoord, dcoord)):
                    if i < len(bootstrap_values) and bootstrap_values[i] > 0:
                        ax.text((ic[1] + ic[2]) / 2, dc[1],
                                f"{bootstrap_values[i]:.0f}%",
                                ha='center', va='bottom', fontsize=8,
                                bbox=dict(facecolor='white', alpha=0.7,
                                          edgecolor='none'))

            plt.tight_layout()
            for fmt, path in [
                ('jpeg', os.path.join(output_dir, "dendrogram.jpg")),
                ('pdf',  os.path.join(output_dir, "dendrogram.pdf")),
            ]:
                plt.savefig(path, format=fmt, dpi=300,
                            bbox_inches='tight', transparent=True)
            plt.close()

            n_spp = len(species_names)
            bio_mat = [
                distance_matrix.values[i, :i + 1].tolist()
                for i in range(n_spp)
            ]
            bio_dm  = DistanceMatrix(names=labels, matrix=bio_mat)
            tree    = DistanceTreeConstructor().upgma(bio_dm)
            Phylo.write(tree,
                        os.path.join(output_dir, "dendrogram.newick"),
                        "newick")

        except PermissionError:
            messagebox.showerror(
                "Error",
                "Permission denied. Close any open files in the output directory.")
        except Exception as dend_exc:
            logging.error("Dendrogram error: %s", dend_exc)
            error_msg = str(dend_exc)
            messagebox.showerror("Dendrogram Error",
                                 f"Failed to generate dendrogram: {error_msg}")

    # ==================== MÉTODOS DE CONFIGURAÇÃO ====================

    def save_primer_settings(self):
        fp = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")])
        if not fp:
            return
        try:
            with open(fp, 'w') as fh:
                json.dump({'min_repeats': self.min_repeats,
                           'primer_params': self.primer_params}, fh, indent=4)
            messagebox.showinfo("Success", "Settings saved!")
            logging.info("Primer settings saved: %s", fp)
        except Exception as save_exc:
            error_msg = str(save_exc)
            messagebox.showerror("Error", f"Failed to save settings: {error_msg}")

    def load_primer_settings(self):
        fp = filedialog.askopenfilename(
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")])
        if not fp:
            return
        try:
            with open(fp) as fh:
                settings = json.load(fh)
            self.min_repeats    = settings.get('min_repeats',    self.min_repeats)
            self.primer_params  = settings.get('primer_params',  self.primer_params)
            self._compile_ssr_patterns()
            self._sync_vars_from_params()
            messagebox.showinfo("Success", "Settings loaded!")
            logging.info("Primer settings loaded: %s", fp)
        except Exception as load_exc:
            error_msg = str(load_exc)
            messagebox.showerror("Error", f"Failed to load settings: {error_msg}")

    def reset_primer_settings(self):
        self.min_repeats = {
            'mono':  {'motif_length': 1, 'min_repeats': 10, 'enabled': True},
            'di':    {'motif_length': 2, 'min_repeats': 5,  'enabled': True},
            'tri':   {'motif_length': 3, 'min_repeats': 4,  'enabled': True},
            'tetra': {'motif_length': 4, 'min_repeats': 3,  'enabled': True},
            'penta': {'motif_length': 5, 'min_repeats': 3,  'enabled': True},
            'hexa':  {'motif_length': 6, 'min_repeats': 3,  'enabled': True},
        }
        self.primer_params = {
            'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 15,
            'PRIMER_MAX_SIZE': 30, 'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 50.0, 'PRIMER_MAX_TM': 70.0,
            'PRIMER_MIN_GC': 30.0, 'PRIMER_MAX_GC': 70.0,
            'PRIMER_OPT_GC_PERCENT': 50.0, 'PRIMER_MAX_POLY_X': 4,
            'PRIMER_MAX_END_GC': 2,
            'PRIMER_PRODUCT_SIZE_RANGE': [[100, 600]],
            'PRIMER_MAX_TEMPLATE_MISPRIMING': 12.0,
            'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING': 24.0,
            'PRIMER_MAX_SELF_ANY': 8.0, 'PRIMER_MAX_SELF_END': 3.0,
            'PRIMER_PAIR_MAX_COMPL_ANY': 8.0,
            'PRIMER_PAIR_MAX_COMPL_END': 3.0,
            'PRIMER_MAX_HAIRPIN_TH': 24.0, 'PRIMER_MAX_DIFF_TM': 2.0,
            'PRIMER_SALT_MONOVALENT': 50.0, 'PRIMER_SALT_DIVALENT': 1.5,
            'PRIMER_DNTP_CONC': 0.8, 'PRIMER_DNA_CONC': 50.0,
        }
        self._compile_ssr_patterns()
        self._sync_vars_from_params()
        messagebox.showinfo("Success", "Settings reset to defaults!")

    def _sync_vars_from_params(self):
        for rtype, info in self.min_repeats.items():
            if hasattr(self, f"{rtype}_var"):
                getattr(self, f"{rtype}_var").set(str(info['min_repeats']))
            if hasattr(self, f"{rtype}_enabled_var"):
                getattr(self, f"{rtype}_enabled_var").set(info['enabled'])
        special = {
            'PRIMER_PRODUCT_SIZE_MIN': self.primer_params[
                'PRIMER_PRODUCT_SIZE_RANGE'][0][0],
            'PRIMER_PRODUCT_SIZE_MAX': self.primer_params[
                'PRIMER_PRODUCT_SIZE_RANGE'][0][1],
        }
        for key, val in {**self.primer_params, **special}.items():
            var_name = f"{key}_var"
            if hasattr(self, var_name):
                getattr(self, var_name).set(str(val))

    def update_parameters(self) -> bool:
        try:
            for rtype in self.min_repeats:
                val = int(getattr(self, f"{rtype}_var").get())
                if val < 1:
                    raise ValueError(f"{rtype} min_repeats must be ≥ 1")
                self.min_repeats[rtype]['min_repeats'] = val
                self.min_repeats[rtype]['enabled'] = \
                    getattr(self, f"{rtype}_enabled_var").get()

            int_params = [
                'PRIMER_MIN_SIZE', 'PRIMER_MAX_SIZE', 'PRIMER_OPT_SIZE',
                'PRIMER_MAX_POLY_X', 'PRIMER_MAX_END_GC',
            ]
            float_params = [
                'PRIMER_MIN_GC', 'PRIMER_MAX_GC', 'PRIMER_OPT_GC_PERCENT',
                'PRIMER_MIN_TM', 'PRIMER_MAX_TM', 'PRIMER_OPT_TM',
                'PRIMER_MAX_DIFF_TM', 'PRIMER_SALT_MONOVALENT',
                'PRIMER_SALT_DIVALENT', 'PRIMER_DNTP_CONC', 'PRIMER_DNA_CONC',
                'PRIMER_MAX_SELF_ANY', 'PRIMER_MAX_SELF_END',
                'PRIMER_PAIR_MAX_COMPL_ANY', 'PRIMER_PAIR_MAX_COMPL_END',
            ]
            for key in int_params:
                val = int(getattr(self, f"{key}_var").get())
                if val <= 0:
                    raise ValueError(f"{key} must be > 0")
                self.primer_params[key] = val

            for key in float_params:
                val = float(getattr(self, f"{key}_var").get())
                if val < 0:
                    raise ValueError(f"{key} must be ≥ 0")
                self.primer_params[key] = val

            ps_min = int(self.PRIMER_PRODUCT_SIZE_MIN_var.get())
            ps_max = int(self.PRIMER_PRODUCT_SIZE_MAX_var.get())
            if ps_min >= ps_max:
                raise ValueError(
                    "Product size MIN must be less than MAX")
            self.primer_params['PRIMER_PRODUCT_SIZE_RANGE'] = [[ps_min, ps_max]]

            if self.primer_params['PRIMER_MIN_SIZE'] >= \
                    self.primer_params['PRIMER_MAX_SIZE']:
                raise ValueError("Primer MIN_SIZE must be < MAX_SIZE")
            if self.primer_params['PRIMER_MIN_TM'] >= \
                    self.primer_params['PRIMER_MAX_TM']:
                raise ValueError("Primer MIN_TM must be < MAX_TM")

            self._compile_ssr_patterns()
            return True

        except ValueError as exc:
            messagebox.showerror(
                "Invalid Parameter",
                f"Please correct the following before running:\n\n{exc}"
            )
            logging.warning("Parameter validation error: %s", exc)
            return False

    # ==================== MÉTODOS DE PROJETO ====================

    def save_project(self):
        fp = filedialog.asksaveasfilename(
            defaultextension=".ssrdev",
            filetypes=[("SSRdev project", "*.ssrdev"),
                       ("All files", "*.*")])
        if not fp:
            return
        try:
            state = {
                'version':       '4.5',
                'project_name':  self.project_entry.get().strip(),
                'fasta_file':    self.file_path.get(),
                'min_repeats':   self.min_repeats,
                'primer_params': self.primer_params,
                'analysis_results': [
                    {k: v for k, v in r.items() if k != 'sequence'}
                    for r in (self.analysis_results or [])
                ],
                'marker_name_map': self.marker_name_map,
                'blast_fasta_file': self.blast_fasta_file,
            }
            with open(fp, 'w') as fh:
                json.dump(state, fh, indent=2)
            messagebox.showinfo("Saved", f"Project saved to:\n{fp}")
            logging.info("Project saved: %s", fp)
        except Exception as save_exc:
            error_msg = str(save_exc)
            messagebox.showerror("Error", f"Could not save project: {error_msg}")

    def load_project(self):
        fp = filedialog.askopenfilename(
            filetypes=[("SSRdev project", "*.ssrdev"),
                       ("All files", "*.*")])
        if not fp:
            return
        try:
            with open(fp) as fh:
                state = json.load(fh)
            self.project_entry.delete(0, tk.END)
            self.project_entry.insert(0, state.get('project_name', ''))
            if state.get('fasta_file'):
                self.file_path.set(state['fasta_file'])
            if state.get('blast_fasta_file'):
                self.blast_fasta_file = state['blast_fasta_file']
                self.blast_fasta_var.set(state['blast_fasta_file'])
                self._load_reference_sequences(state['blast_fasta_file'])
            self.min_repeats   = state.get('min_repeats',   self.min_repeats)
            self.primer_params = state.get('primer_params', self.primer_params)
            self.marker_name_map = state.get('marker_name_map', {})
            self._compile_ssr_patterns()
            self._sync_vars_from_params()
            if state.get('analysis_results'):
                self.analysis_results = state['analysis_results']
                self.display_results()
            messagebox.showinfo("Loaded", f"Project loaded from:\n{fp}")
            logging.info("Project loaded: %s", fp)
        except Exception as load_exc:
            error_msg = str(load_exc)
            messagebox.showerror("Error", f"Could not load project: {error_msg}")

    # ==================== MÉTODOS AUXILIARES ====================

    def browse_file(self):
        fp = filedialog.askopenfilename(
            filetypes=[('FASTA files', '*.fasta *.fa'),
                       ('All files', '*.*')])
        if fp:
            self.file_path.set(fp)
            self.input_directory = os.path.dirname(fp)

    def browse_fasta(self):
        fp = filedialog.askopenfilename(
            filetypes=[("FASTA files", "*.fasta *.fa *.fna"),
                       ("All files", "*.*")])
        if fp:
            self.fasta_file.set(fp)
            base = os.path.splitext(os.path.basename(fp))[0]
            self.output_file.set(
                os.path.join(os.path.dirname(fp),
                             f"{base}_validation_results.xlsx"))

    def browse_excel(self):
        fp = filedialog.askopenfilename(
            filetypes=[("Excel files", "*.xlsx *.xls"),
                       ("All files", "*.*")])
        if fp:
            self.excel_file.set(fp)

    def get_sequences(self) -> list:
        sequences = []
        fp = self.file_path.get()
        if fp:
            try:
                sequences.extend(list(SeqIO.parse(fp, 'fasta')))
            except Exception as exc:
                error_msg = str(exc)
                messagebox.showerror("Error", f"Error reading file:\n{error_msg}")
                return []
        text = self.seq_text.get('1.0', tk.END).strip()
        if text:
            try:
                sequences.extend(list(SeqIO.parse(StringIO(text), 'fasta')))
            except Exception as exc:
                error_msg = str(exc)
                messagebox.showerror("Error", f"Error parsing text:\n{error_msg}")
                return []
        if not sequences:
            messagebox.showwarning("Warning", "No sequences provided")
        return sequences

    def show_page(self, page_name: str):
        pages = {"Prospecting": 0, "Primer Settings": 1,
                 "Characterization": 2, "Validation": 3}
        self.notebook.select(pages[page_name])

    def clear_results(self):
        self.analysis_results = None
        for w in self.results_frame.winfo_children():
            w.destroy()
        self.progress['value'] = 0
        self.status.config(text="Ready")

    def select_family(self):
        families = {
            "Plants": [
                "Acanthaceae","Aceraceae","Acoraceae","Adoxaceae","Aizoaceae",
                "Alismataceae","Amaranthaceae","Amaryllidaceae","Anacardiaceae",
                "Apiaceae","Apocynaceae","Araceae","Araliaceae","Arecaceae",
                "Asparagaceae","Asteraceae","Balsaminaceae","Berberidaceae",
                "Betulaceae","Bignoniaceae","Boraginaceae","Brassicaceae",
                "Bromeliaceae","Burseraceae","Cactaceae","Campanulaceae",
                "Cannabaceae","Caprifoliaceae","Caryophyllaceae","Celastraceae",
                "Chenopodiaceae","Clusiaceae","Colchicaceae","Commelinaceae",
                "Convolvulaceae","Cornaceae","Crassulaceae","Cucurbitaceae",
                "Cupressaceae","Cyperaceae","Dipsacaceae","Droseraceae",
                "Ebenaceae","Elaeagnaceae","Ericaceae","Euphorbiaceae",
                "Fabaceae","Fagaceae","Gentianaceae","Geraniaceae",
                "Gesneriaceae","Grossulariaceae","Hamamelidaceae",
                "Hydrangeaceae","Hypericaceae","Iridaceae","Juglandaceae",
                "Lamiaceae","Lauraceae","Lentibulariaceae","Liliaceae",
                "Linaceae","Lythraceae","Magnoliaceae","Malvaceae",
                "Melastomataceae","Meliaceae","Menispermaceae","Moraceae",
                "Myrtaceae","Nyctaginaceae","Nymphaeaceae","Oleaceae",
                "Onagraceae","Orchidaceae","Orobanchaceae","Oxalidaceae",
                "Paeoniaceae","Papaveraceae","Passifloraceae","Phyllanthaceae",
                "Pinaceae","Plantaginaceae","Platanaceae","Plumbaginaceae",
                "Poaceae","Polemoniaceae","Polygalaceae","Polygonaceae",
                "Portulacaceae","Primulaceae","Proteaceae","Ranunculaceae",
                "Rhamnaceae","Rosaceae","Rubiaceae","Rutaceae","Salicaceae",
                "Santalaceae","Sapindaceae","Sapotaceae","Saxifragaceae",
                "Scrophulariaceae","Simaroubaceae","Solanaceae","Tamaricaceae",
                "Theaceae","Thymelaeaceae","Tiliaceae","Tropaeolaceae",
                "Ulmaceae","Urticaceae","Verbenaceae","Violaceae","Vitaceae",
                "Zingiberaceae","Zygophyllaceae",
            ],
            "Animals": [
                "Acanthocephala","Acoelomorpha","Annelida","Arthropoda",
                "Brachiopoda","Bryozoa","Chaetognatha","Chordata","Cnidaria",
                "Ctenophora","Cycliophora","Echinodermata","Entoprocta",
                "Gastrotricha","Gnathostomulida","Hemichordata","Kinorhyncha",
                "Loricifera","Micrognathozoa","Mollusca","Nematoda",
                "Nematomorpha","Nemertea","Onychophora","Orthonectida",
                "Phoronida","Placozoa","Platyhelminthes","Porifera",
                "Priapulida","Rhombozoa","Rotifera","Tardigrada",
                "Xenacoelomorpha",
            ],
            "Fungi": [
                "Ascomycota","Basidiomycota","Blastocladiomycota",
                "Chytridiomycota","Cryptomycota","Glomeromycota",
                "Microsporidia","Neocallimastigomycota","Zygomycota",
                "Agaricaceae","Amanitaceae","Boletaceae","Cantharellaceae",
                "Clavariaceae","Cortinariaceae","Entolomataceae",
                "Ganodermataceae","Hydnaceae","Hygrophoraceae","Inocybaceae",
                "Marasmiaceae","Meripilaceae","Mycenaceae","Omphalotaceae",
                "Physalacriaceae","Pleurotaceae","Pluteaceae","Polyporaceae",
                "Psathyrellaceae","Russulaceae","Schizophyllaceae",
                "Sclerodermataceae","Strophariaceae","Tricholomataceae",
                "Tuberaceae",
            ],
            "Protists": [
                "Amoebozoa","Apicomplexa","Choanoflagellata","Ciliophora",
                "Dinoflagellata","Diplomonadida","Euglenozoa","Foraminifera",
                "Heterolobosea","Kinetoplastida","Mycetozoa","Parabasalia",
                "Radiolaria","Stramenopiles","Alveolata","Cercozoa",
                "Cryptophyta","Glaucophyta","Haptophyta","Rhizaria",
            ],
            "Algae": [
                "Chlorophyta","Charophyta","Rhodophyta","Phaeophyta",
                "Bacillariophyta","Chrysophyta","Dinophyta","Euglenophyta",
                "Xanthophyta","Prasinophyta","Trebouxiophyta","Ulvophyta",
                "Zygnematophyta",
            ],
            "Bacteria": [
                "Acidobacteria","Actinobacteria","Aquificae","Bacteroidetes",
                "Chlamydiae","Chlorobi","Chloroflexi","Cyanobacteria",
                "Deinococcus-Thermus","Firmicutes","Fusobacteria",
                "Gemmatimonadetes","Nitrospirae","Planctomycetes",
                "Proteobacteria","Spirochaetes","Thermotogae",
                "Verrucomicrobia",
            ],
            "Archaea": [
                "Crenarchaeota","Euryarchaeota","Korarchaeota",
                "Nanoarchaeota","Thaumarchaeota",
            ],
        }

        win = tk.Toplevel(self.root)
        win.title("Select Family")
        win.geometry("600x700")

        mf = ttk.Frame(win, padding=10)
        mf.pack(fill=tk.BOTH, expand=True)

        sf  = ttk.Frame(mf)
        sf.pack(fill=tk.X, pady=(0, 10))
        ttk.Label(sf, text="Search:").pack(side=tk.LEFT, padx=(0, 5))
        search_var = tk.StringVar()
        se = ttk.Entry(sf, textvariable=search_var)
        se.pack(side=tk.LEFT, fill=tk.X, expand=True)
        se.focus()

        tf  = ttk.Frame(mf)
        tf.pack(fill=tk.BOTH, expand=True)
        tree = ttk.Treeview(tf, columns=('Family',),
                            show='tree headings', height=20)
        tree.heading('#0', text='Group')
        tree.heading('Family', text='Family')
        tree.column('#0', width=150, minwidth=150)
        tree.column('Family', width=300, minwidth=300)

        vsb = ttk.Scrollbar(tf, orient="vertical", command=tree.yview)
        hsb = ttk.Scrollbar(tf, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        tree.grid(row=0, column=0, sticky='nsew')
        vsb.grid(row=0, column=1, sticky='ns')
        hsb.grid(row=1, column=0, sticky='ew')
        tf.grid_rowconfigure(0, weight=1)
        tf.grid_columnconfigure(0, weight=1)

        for group, flist in families.items():
            gid = tree.insert('', 'end', text=group, values=('',))
            for fam in sorted(flist):
                tree.insert(gid, 'end', text='', values=(fam,))

        def _filter(*_):
            term = search_var.get().lower()
            for item in tree.get_children():
                tree.item(item, open=False)
                for child in tree.get_children(item):
                    if term in tree.item(child, 'values')[0].lower():
                        tree.item(item, open=True)
                        tree.see(child)
                        break
        search_var.trace_add('write', _filter)

        bf = ttk.Frame(mf)
        bf.pack(fill=tk.X, pady=(10, 0))

        def _select():
            sel = tree.selection()
            if sel:
                item = tree.item(sel[0])
                if item['values'][0]:
                    parent = tree.parent(sel[0])
                    grp    = tree.item(parent, 'text')
                    fam    = item['values'][0]
                    self.family_var.set(f"{grp}: {fam}")
                    win.destroy()
                else:
                    messagebox.showwarning(
                        "Warning", "Please select a family, not a group")

        ttk.Button(bf, text="Select", command=_select,
                   style="Accent.TButton").pack(side=tk.RIGHT, padx=5)
        ttk.Button(bf, text="Cancel",
                   command=win.destroy).pack(side=tk.RIGHT, padx=5)

        win.transient(self.root)
        win.grab_set()
        win.wait_window()


if __name__ == "__main__":
    root = tk.Tk()
    app  = SSRdevApp(root)
    root.mainloop()
