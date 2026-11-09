#!/usr/bin/env python3
"""Generate docs/RNA.html and docs/DNA.html (self-contained) from the content below.
Version status is read from the repo tree (archive/<name>/ file count)."""
import html, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
E = html.escape

def nver(name):
    stem = name.rsplit(".", 1)[0]
    d = os.path.join(ROOT, "archive", stem)
    n = len(os.listdir(d)) if os.path.isdir(d) else 0
    return n

CSS = """body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:0;padding:0 14px 40px;max-width:900px;margin:auto;line-height:1.45;color:#222}
h1{font-size:1.5em}h2{border-bottom:2px solid #4a6fa5;padding-top:12px;font-size:1.25em}h3{font-size:1.05em;margin-bottom:4px}
code,pre{font-family:Menlo,Consolas,monospace;font-size:.85em}pre{background:#f3f5f8;padding:8px;overflow-x:auto;border-radius:4px}
.step{border:1px solid #d5dbe5;border-radius:6px;padding:8px 12px;margin:10px 0}.step p{margin:4px 0}
.k{font-weight:600;color:#4a6fa5}.gloss{background:#fff8e1;padding:2px 6px;border-radius:3px}
.ver{font-size:.85em;color:#555}nav a{margin-right:8px;white-space:nowrap}svg{max-width:100%;height:auto}
.unclear{color:#b3261e;font-weight:600}"""

def svg(stages):
    w = 150; h = 46; gap = 18
    per_row = 3
    rows = (len(stages) + per_row - 1) // per_row
    W = per_row * (w + gap) + gap; H = rows * (h + 30) + 10
    out = [f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" font-family="sans-serif" font-size="11">',
           '<defs><marker id="a" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#4a6fa5"/></marker></defs>']
    pos = []
    for i, s in enumerate(stages):
        r, c = divmod(i, per_row)
        if r % 2: c = per_row - 1 - c
        x = gap + c * (w + gap); y = 10 + r * (h + 30); pos.append((x, y))
        out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" fill="#e8eef8" stroke="#4a6fa5"/>')
        lines = s.split("|")
        for j, t in enumerate(lines):
            out.append(f'<text x="{x+w/2}" y="{y+18+14*j}" text-anchor="middle">{E(t)}</text>')
    for i in range(len(stages) - 1):
        (x1, y1), (x2, y2) = pos[i], pos[i + 1]
        if y1 == y2:
            a = (x1 + w, y1 + h / 2, x2, y2 + h / 2) if x2 > x1 else (x1, y1 + h / 2, x2 + w, y2 + h / 2)
        else:
            a = (x1 + w / 2, y1 + h, x2 + w / 2, y2)
        out.append(f'<line x1="{a[0]}" y1="{a[1]}" x2="{a[2]}" y2="{a[3]}" stroke="#4a6fa5" stroke-width="1.5" marker-end="url(#a)"/>')
    out.append("</svg>")
    return "\n".join(out)

def page(title, fam, intro, stages, sections, notes):
    o = [f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{E(title)}</title><style>{CSS}</style></head><body>']
    o.append(f"<h1>{E(title)}</h1>{intro}")
    o.append('<p><b>Data:</b> not in this repo. See <code>FILES.html</code> in the Yunshi_data folder on the USB drive.</p>')
    o.append("<nav>" + "".join(f'<a href="#{i}">{E(t)}</a>' for i, t, _ in sections) + "</nav>")
    o.append("<h2>Workflow at a glance</h2>" + svg(stages))
    o.append("<p>Paths below are relative to the repo; scripts are run inside the simulation or analysis directory, not from the repo root. Hard-coded cluster paths are listed in <a href='PORTABILITY.md'>PORTABILITY.md</a>. Versions: <a href='VERSIONS.md'>VERSIONS.md</a>.</p>")
    for i, t, steps in sections:
        o.append(f'<h2 id="{i}">{E(t)}</h2>')
        for st in steps:
            path, when, io, cmd, gloss = st
            name = os.path.basename(path.split()[0])
            n = nver(name) if "." in name and not path.startswith("templates") else 0
            ver = f"canonical = newest copy; {n} older distinct version(s) in <code>archive/{E(name.rsplit('.',1)[0])}/</code>" if n else "single version"
            o.append(f'<div class="step"><h3><code>{E(path)}</code></h3>')
            o.append(f'<p><span class="k">When:</span> {when}</p>')
            if gloss: o.append(f'<p class="gloss">{gloss}</p>')
            o.append(f'<p><span class="k">In / out:</span> {io}</p>')
            if cmd: o.append(f"<pre>{E(cmd)}</pre>")
            o.append(f'<p class="ver">Version: {ver}.</p></div>')
    o.append("<h2>Unclear / caveats</h2><ul>" + "".join(f"<li>{n}</li>" for n in notes) + "</ul>")
    o.append("</body></html>")
    return "\n".join(o)

# ---------------- RNA ----------------
RNA_SECTIONS = [
 ("setup", "1. Setup of a run", [
  ("rna/setup/setup_run.py", "First step for a new structure. Run in a directory holding exactly one <code>.top</code> and one <code>.dat/.oxdna</code>.",
   "In: topology + start configuration. Out: <code>&lt;prefix&gt;_relax/_md/_analysis</code> dirs, copies of <code>input_pre_relax</code>, <code>input_cuda_relax</code>, <code>plot_energy.ipynb</code>, and <code>cmds_relax</code>.",
   "python setup_run.py\n# cmds_relax then contains:\n# /home/dfan/oxDNA/build/bin/oxDNA input_pre_relax\n# sleep 5\n# /home/dfan/oxDNA/build/bin/oxDNA input_cuda_relax",
   "oxRNA: coarse-grained RNA model, one bead-with-sites per nucleotide."),
  ("rna/setup/replacer.py", "After setup_run: fill in file names in the input files.", "In: input* files. Out: <code>replace_me.dat/.top</code> replaced by the real names.", "python replacer.py", ""),
  ("rna/setup/special_replacer.py", "Variant of replacer for non-standard placeholder cases. <span class='unclear'>Unclear</span> exact use from code alone.", "In/out: input files edited in place.", "", ""),
  ("rna/setup/mt_creater.py", "To link monomers into an N-mer (8mer): writes mutual_trap forces.", "Out: <code>forces.txt</code> with stiff springs, r0 = 1.2, PBC=1. MONOMER_LENGTH = 244.", "python mt_creater.py   # edit k-mer settings in the script", "Mutual trap: harmonic spring between two nucleotides."),
  ("rna/setup/al_mt_creater.py", "Alternative (aligned) mutual-trap generator.", "Out: forces file.", "", ""),
  ("rna/setup/multimerforce.py", "Older/alternate multimer force generator.", "Out: forces file.", "", ""),
  ("rna/setup/force_remover.py", "Remove force lines (e.g. before a production run without the linkers).", "In: forces file. Out: edited forces file.", "", ""),
  ("rna/setup/ntremover.py", "Remove nucleotides from a topology/configuration pair. <span class='unclear'>Unclear</span> argument order; read the script header.", "In: .top/.dat.", "", ""),
  ("templates/rna/input_pre_relax, input_cuda_relax, input_md", "Templates the setup scripts copy; edit temperature, salt, steps.", "Used by setup_run.py.", "", "")]),
 ("relax", "2. MC / MD relaxation", [
  ("templates/rna/input_pre_relax", "Monte Carlo pre-relaxation on CPU, to remove overlaps.", "Out: relaxed configuration for the next stage.", "/home/dfan/oxDNA/build/bin/oxDNA input_pre_relax", "MC: random moves accepted by Metropolis criterion."),
  ("templates/rna/input_cuda_relax", "Short MD relaxation on GPU (CUDA).", "In: MC output. Out: relaxed conf + energy file (plot with notebooks/rna/plot_energy.ipynb).", "/home/dfan/oxDNA/build/bin/oxDNA input_cuda_relax", "MD: integrates Newton equations with a thermostat."),
  ("jobs/rna/queue.sh", "Submit relaxation commands on Slurm as an array (one line of cmds per task).", "In: file <code>cmds</code>.", "sbatch queue.sh", "")]),
 ("prod", "3. Production MD", [
  ("templates/rna/input_md", "Production MD input template (set steps, print intervals, T).", "Out: trajectory.dat, energy.dat, last_conf.dat.", "oxDNA input_md", ""),
  ("jobs/rna/queue_25_prod.sh", "Submit production runs on Slurm (GPU partition).", "In: cmds file.", "sbatch queue_25_prod.sh", "")]),
 ("analysis", "4. Analysis of trajectories", [
  ("rna/setup/setup_analysis.py", "Prepare the analysis directory; writes the standard command list.", "Out: <code>cmds.sh</code>.", "python setup_analysis.py\n# cmds.sh:\noat mean -p 4 -o mean.dat -d devs.json *trj.dat\noat align -p 4 *trj.dat align.dat\npython hb_cal.py *_input.txt *trj.dat kissingloop.csv --output KL_HB.csv\npython total_hbond_cal.py *input.txt *trj.dat --output Total_HB.csv\npython signed_helical_bent_twist_cal.py align.dat helixces.txt", "oat: oxDNA_analysis_tools, command-line toolkit (mean structure, alignment, RMSF)."),
  ("rna/setup/prepare_oat_workflow.py", "Alternative that builds oat bond analysis from the forces.", "Out: cmds.sh.", "oat generate_forces input_pre_relax last_conf.dat\noat forces2pairs forces.txt\noat bond_analysis *input.txt *_trj.dat pairs.txt -p 4", ""),
  ("rna/analysis/hb_cal.py", "Count hydrogen bonds among a chosen set of nucleotides (e.g. kissing loop).", "In: input, trajectory, index csv. Out: csv. Occupancy % = bonds / (n/2) x 100.", "python hb_cal.py input.txt trj.dat kissingloop.csv -p 4 --output KL_HB.csv", "Hydrogen bond (hb_list observable): base-pair contact energy below a cutoff."),
  ("rna/analysis/total_hbond_cal.py", "Same for all nucleotides.", "Out: Total_HB.csv.", "python total_hbond_cal.py input.txt trj.dat --output Total_HB.csv", ""),
  ("rna/analysis/calculate_bond_instability.py", "Per-bond stability from a bond-occupancy table.", "Out: Bond Instability Index 4p(1-p) per bond.", "", "BII: 0 for always/never bonded, 1 for bonded half the time."),
  ("rna/analysis/log_transformer.py", "Rescale BII for plotting.", "log(BII+1e-6) then min-max to [0,1].", "", ""),
  ("rna/analysis/plot_bond_instability.py", "Plot BII.", "", "", ""),
  ("rna/analysis/calculate_bond.py", "Bond summary helper. <span class='unclear'>Unclear</span> how it differs from hb_cal.py.", "", "", ""),
  ("rna/analysis/signed_helical_bent_twist_cal.py", "Main bend/twist measurement per frame (current).", "In: aligned trajectory + junctions file (1-based inclusive ranges, 3 or 5 fields). Out: <code>twist_bend.csv</code> (t, twist_deg, bend_deg).", "python signed_helical_bent_twist_cal.py align.dat helixces.txt", "Segments LH1, SH1, SH2, LH4 (long/short helices). Bend: angle between SVD axes of LH1 and LH4. Twist: angle between radial vectors of SH1 and SH2 about the long-helix axis, signed right-handed. Bend sign from a reference normal fixed on frame 1."),
  ("rna/analysis/helical_bent_twist_cal.py", "Unsigned predecessor.", "Out: twist/bend csv.", "", ""),
  ("rna/analysis/helix_bt.py, twist_bent.py, tb2.py", "Earlier bend/twist variants, kept for older runs.", "", "", ""),
  ("rna/analysis/distance_cal.py", "Distance between two nucleotides over time.", "Out: distance file.", "python distance_cal.py trj.dat nt1 nt2", ""),
  ("rna/analysis/dis_bent_plotter.py, dis_twist_plotter.py, dis_bond_plotter.py, dis_energy_plotter.py", "Plot distance vs bend / twist / bonds / energy.", "In: distance file + the other file.", "python dis_bent_plotter.py distance_file bend_file", ""),
  ("rna/analysis/plot_energydat.py", "Plot energy.dat.", "", "", ""),
  ("rna/analysis/persistence_length.py", "Persistence length of a helical segment from tangent correlations.", "In: trajectory, input, csv of segment. Out: Lp (nm) + plot.", "python persistence_length.py traj.dat input.txt segment.csv -p 4", "Persistence length Lp: distance over which direction memory decays; fit ln C(n) = -n/Lp on first third of offsets, x 0.8518 nm per sim unit."),
  ("rna/multimer/multimer_tb.py", "Bend/twist per monomer of a k-mer.", "Out: monomer_i.csv, average.csv.", "python multimer_tb.py traj.dat junctions.txt 8", ""),
  ("rna/multimer/pairwise_tb.py", "Pairwise version relative to mean reference.", "In: mean structure file.", "python pairwise_tb.py traj.dat junctions.txt 8 average.dat", "")]),
 ("umbrella", "5. Umbrella sampling and WHAM", [
  ("rna/umbrella/setup_umb.py", "Set up umbrella windows along an end-to-end reaction coordinate.", "In: parameters.txt (topology, xi_start, xi_end, n_windows, stiff, anchor_a, anchor_b, steps, template, T, salt). Out: <code>window_NN/equil|prod</code> with input, forces, check_equil.py, plus run_all.sh.", "python setup_umb.py parameters.txt\n# e.g. xi 7.5-10.2, 20 windows, stiff 100, anchors 16/226", "Umbrella sampling: each window holds the coordinate near a target with a harmonic bias, so rare states are sampled."),
  ("rna/umbrella/setup_prod.py", "Regenerate production-only scripts.", "Out: run_prod.sh.", "", ""),
  ("rna/umbrella/tune_stiff.py", "Change stiffness / steps / print intervals in production inputs.", "Out: run_production_only.sh.", "python tune_stiff.py NEW_STIFF NEW_STEPS NEW_PRINT_CONF NEW_PRINT_ENERGY [OXDNA_BIN]", ""),
  ("rna/umbrella/run_all.sh, run_prod.sh, run_production_only.sh", "Auto-generated drivers: chained equilibrations then production. Contain absolute cluster paths.", "", "bash run_all.sh", ""),
  ("rna/umbrella/check_equil.py", "Per window: check equilibration of the coordinate. Auto-generated per window (only ANCHOR_A, ANCHOR_B, XI0 differ).", "", "", ""),
  ("rna/umbrella/check_dis.py", "Quick distance check of a configuration.", "", "python check_dis.py conf.dat anchor_a anchor_b", ""),
  ("rna/umbrella/wham_cal.py, wham_com.py, wham_calculate.py", "Combine windows into a free-energy profile.", "Out: PMF (.dat) and plot.", "python wham_cal.py   # run in umbrella dir", "WHAM: iterative reweighting of biased histograms into one unbiased PMF. Bias 0.5 k (x-r0)^2, beta=1, 80 bins, converge max|df|<1e-6 (<=10000 iters), F = -ln P, min shifted to 0, x 0.8518 nm. wham_com uses centre-of-mass distance."),
  ("rna/umbrella/force_from_pmf.py", "Force profile from PMF.", "f = -dF/dx x 4.1 pN per kBT/nm.", "", ""),
  ("rna/umbrella/compute_barrier_forces.py", "Barrier height, closed/open minima, DeltaF, mean downhill forces.", "", "", ""),
  ("rna/umbrella/2state.py, state_average.py", "Two-state analysis: open = 0 HB, closed = 6 HB.", "DeltaG = -kBT ln(P_closed/P_open).", "", ""),
  ("rna/umbrella/plot_window.py, pwindow.py, plot_complete.py, plot_completed.py, plot_com_by_window.py", "Plots per window and overall.", "", "", "")]),
 ("ffs", "6. Forward flux sampling (FFS)", [
  ("rna/ffs/setup_kl_ffs.py", "Set up an FFS run for a kissing loop.", "Copies op.txt, creates FLUX and N_shoot dirs, writes cmd_ffs.", "python setup_kl_ffs.py X", "FFS: estimates rates of rare transitions by shooting trajectories between successive order-parameter interfaces."),
  ("rna/ffs/ffs_flux.py", "Stage 1: flux through the first interface.", "Template: templates/rna/ffs/input (sim_type=FFS_MD), op.txt.", "python ffs_flux.py -n 100 -s 1 -c 8", ""),
  ("rna/ffs/ffs_shoot.py", "Stage 2: shoot from one interface to the next.", "", "python ffs_shoot.py -n 100 -s 1 -c 8 -k 0", ""),
  ("rna/ffs/extract_conf_hb.py", "Extract configurations / HB counts from FFS output. <span class='unclear'>Unclear</span> details.", "", "", "")]),
 ("ramp", "7. Thermal ramp", [
  ("rna/thermal_ramp/thermal_ramp_setup.py", "Make per-step inputs for a temperature ramp.", "Out: inputMD_stepN_XXc.", "python thermal_ramp_setup.py -s 37 -e 4 -t 10 -input input_md", ""),
  ("rna/thermal_ramp/run_thermal_ramp.sh", "Run steps 1..10 (37 C down to 4 C) in sequence.", "", "bash run_thermal_ramp.sh", ""),
  ("rna/thermal_ramp/wrap_thermal_ramp.py", "Join step trajectories and align.", "Runs <code>oat align</code>.", "", "")]),
 ("plot", "8. Plotting notebooks", [
  ("notebooks/rna/plot_energy.ipynb", "Energy of MC and MD relaxation.", "In: mc_relax_energy.dat, md_relax_energy.dat.", "", ""),
  ("notebooks/rna/plot_analysis_results.ipynb", "Bend/twist time series, histograms, HB.", "In: twist_bend.csv.", "", ""),
  ("notebooks/rna/plot_multimer.ipynb", "Per-monomer twist/bend (k-mer).", "In: monomer_i.csv, average.csv.", "", ""),
  ("notebooks/rna/cross_plotter.ipynb, plot_monos.ipynb, figure_template.ipynb", "Cross-design KDE comparisons and figure styling. figure_template has a hard-coded local Downloads path.", "", "", "")]),
]
RNA_STAGES = ["Setup|setup_run.py","MC relax|oxDNA CPU","MD relax|oxDNA CUDA","Production MD|input_md","Analysis|oat, hb_cal, bend/twist","Umbrella+WHAM|setup_umb, wham","FFS / thermal ramp|rates, T-dependence","Plots|notebooks"]
RNA_NOTES = ["Which of the multimer-force scripts (mt_creater, al_mt_creater, multimerforce) was used for a given run is <span class='unclear'>unclear</span> from code alone.",
 "umbrella K1IVT8mer and us_k1 directories were filed under RNA by name; confirm.",
 "Exact numerical thresholds in hb_cal come from the oxDNA hb_list observable and the input file, not from the script."]

# ---------------- DNA ----------------
DNA_SECTIONS = [
 ("setup", "1. Setup of a run", [
  ("dna/setup/setup_run.py", "First step for a new nanotweezer structure; same logic as RNA version.", "Out: relax/md/analysis dirs + cmds_relax. Three versions differ only in absolute paths (/home/donglifan vs /home/dfan).", "python setup_run.py", "oxDNA: coarse-grained DNA model."),
  ("dna/setup/replacer.py", "Fill file names into input files.", "", "python replacer.py", ""),
  ("dna/setup/setup_equil_cmds.py", "Command list for equilibration stages.", "", "", ""),
  ("dna/setup/multimer_forces.py", "Linker forces to build a multimer (MONOMER_LENGTH 258; A/B/C/D pair positions).", "Out: forces file.", "", "Mutual trap: harmonic spring between two nucleotides."),
  ("dna/setup/setup_analysis.py, prepare_oat_workflow.py", "Write analysis command lists (as for RNA).", "Out: cmds.sh.", "", "")]),
 ("relax", "2. MC / MD relaxation", [
  ("templates/dna/input_pre_relax", "MC pre-relax on CPU.", "", "oxDNA input_pre_relax", "MC: Metropolis random moves."),
  ("templates/dna/input_cuda_relax", "CUDA MD relax.", "", "oxDNA input_cuda_relax", ""),
  ("jobs/dna/mcjob.sh, mdjob.sh, cpu_mono.sh, queue_relax.sh, equil_job.sh", "Slurm jobs for these stages.", "", "sbatch mdjob.sh", "")]),
 ("prod", "3. Production MD", [
  ("templates/dna/input_md", "Production template.", "", "", ""),
  ("jobs/dna/prod_job.sh, queue_25_prod.sh", "Slurm production jobs.", "", "sbatch prod_job.sh", "")]),
 ("analysis", "4. Analysis", [
  ("dna/analysis/hb_cal.py, total_hbond_cal.py", "Hydrogen-bond counts for a nucleotide set / all.", "Same CLI as RNA.", "python hb_cal.py input.txt trj.dat idx.csv --output HB.csv", "Hydrogen bond (hb_list observable)."),
  ("dna/analysis/signed_helical_bent_twist_cal.py", "Signed bend and twist (current).", "Out: twist_bend.csv.", "python signed_helical_bent_twist_cal.py align.dat junctions.txt", "Bend: angle between long-helix axes; twist about the long-helix axis."),
  ("dna/analysis/angle_cal.py", "Angle between two helices by PCA axes.", "In: trj.dat helices.csv. Out: angles.csv.", "python angle_cal.py trj.dat helices.csv", "PCA axis: principal direction of the nucleotide positions."),
  ("dna/analysis/distance_cal.py, distance_two_sets.py", "Distance between two nucleotides / COM of two groups.", "", "python distance_two_sets.py a.csv b.csv x.dat", "COM: centre of mass."),
  ("dna/analysis/helical_bent_twist_cal.py, helix_bt.py, twist_bent.py, tb2.py", "Earlier bend/twist variants.", "", "", ""),
  ("dna/analysis/calculate_bond.py, calculate_bond_instability.py", "Bond occupancy and Bond Instability Index.", "", "", "BII = 4p(1-p)."),
  ("dna/analysis/dis_*_plotter.py, plot_energydat.py", "Plots vs distance / energy.", "", "", ""),
  ("dna/multimer/multimer_tb.py", "Per-monomer bend/twist for k-mers.", "", "python multimer_tb.py traj.dat junctions.txt k", "")]),
 ("umbrella", "5. Umbrella sampling and WHAM", [
  ("dna/umbrella/setup_umb.py", "Windows along a COM-COM or nucleotide-pair coordinate.", "Out: window_NN/equil|prod, a_coord.txt, b_coord.txt (COM groups).", "python setup_umb.py parameters.txt", "Umbrella sampling: harmonic bias windows."),
  ("dna/umbrella/setup_umbrella.py, setup_prod.py", "Alternative setup / production regeneration.", "", "", ""),
  ("dna/umbrella/tune_stiff.py", "Rewrite stiffness/steps.", "", "python tune_stiff.py NEW_STIFF NEW_STEPS NEW_PRINT_CONF NEW_PRINT_ENERGY", ""),
  ("dna/umbrella/run_all.sh, run_production_only.sh, try.sh", "Generated drivers. run_all.sh holds /Users/a000000/db paths from the umbrella_DNA_opening version.", "", "bash run_all.sh", ""),
  ("dna/umbrella/check_equil.py, check_dis.py", "Equilibration and distance checks (auto-generated per window).", "", "", ""),
  ("dna/umbrella/wham_cal.py, wham_com.py, wham_calculate.py", "WHAM; wham_calculate uses COM from a_coord/b_coord (80 bins, max 10000 frames per window, 0.8518 nm per sim unit).", "Out: PMF.", "python wham_calculate.py", "WHAM: reweights biased windows into one PMF."),
  ("dna/umbrella/tb_wham.py", "Apparent PMFs of twist and bend over windows (no reweighting).", "In: helices.csv. Out: pmf_twist.dat, pmf_bend.dat + violin plots.", "python tb_wham.py helices.csv", "T = 300 K, 80 bins."),
  ("dna/umbrella/force_from_pmf.py, compute_barrier_forces.py", "Force profile and barrier analysis.", "", "", ""),
  ("dna/umbrella/2state.py, state_average.py, delta_energy.py", "Two-state populations and DeltaG.", "", "", "DeltaG = -kBT ln(P_closed/P_open)."),
  ("dna/umbrella/plot_*.py, pwindow.py", "Plots.", "", "", "")]),
 ("ramp", "6. Thermal ramp", [
  ("dna/thermal_ramp/thermal_ramp_setup.py, wrap_thermal_ramp.py", "Per-step inputs and joined trajectory (same as RNA).", "", "python thermal_ramp_setup.py -s 37 -e 4 -t 10 -input input_md", "")]),
 ("plot", "7. Plotting notebooks", [
  ("notebooks/dna/plot_analysis_results.ipynb, plot_energy.ipynb, plot_curves.ipynb", "Results, energy and angle distribution plots.", "In: angles.csv, twist_bend.csv.", "", "")]),
]
DNA_STAGES = ["Setup|setup_run.py","MC relax|oxDNA CPU","MD relax|oxDNA CUDA","Production MD|input_md","Analysis|HB, angle, bend/twist","Umbrella+WHAM|COM coordinate","Thermal ramp|T-dependence","Plots|notebooks"]
DNA_NOTES = ["No FFS scripts exist for DNA in the sources.","dna/umbrella/try.sh purpose is <span class='unclear'>unclear</span>.","umbrella_DNA/us_k1 and umbrella_k1 were classified as RNA by name; confirm."]

def main():
    d = os.path.join(ROOT, "docs"); os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "RNA.html"), "w").write(page("RNA kissing-loop simulations: scripts by workflow", "rna", "<p>oxRNA simulations of kissing-loop designs 2H, K1, T3, IVT (monomer and 8mer). Follow sections top to bottom.</p>", RNA_STAGES, RNA_SECTIONS, RNA_NOTES))
    open(os.path.join(d, "DNA.html"), "w").write(page("DNA nanotweezer simulations: scripts by workflow", "dna", "<p>oxDNA simulations of nanotweezers (nanotw, T3IVT) with umbrella sampling.</p>", DNA_STAGES, DNA_SECTIONS, DNA_NOTES))
main()
