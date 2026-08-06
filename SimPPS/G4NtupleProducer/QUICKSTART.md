# G4NtupleProducer Quick Reference

## Setup (One-time)

```bash
cd ~/PPS/CMSSW_14_0_0
eval "$(scramv1 runtime -sh)"
scramv1 b SimPPS/G4NtupleProducer -j 4
```

## 1. Generate G4 Simulation

```bash
cd ~/PPS/pps2_simulation_prototype

# Option A: Build and run (first time)
nix develop --command bash -c "
    cd build && cmake .. && make -j4 && cd ..
    ./build/PPS2sim macros/run.mac
"

# Option B: Just run (if already built)
nix develop --command bash -c "./build/PPS2sim macros/run.mac"
```

Output: `test_output.root` (check your macro's `/output/filename` directive)

## 2. Run CMSSW Reconstruction

```bash
cd ~/PPS/CMSSW_14_0_0

# Basic (use defaults: arm=0, station=0, rp=3)
cmsRun src/SimPPS/G4NtupleProducer/test/template_cfg.py \
    inputFile=../pps2_simulation_prototype/test_output.root \
    maxEvents=-1

# With custom parameters
cmsRun src/SimPPS/G4NtupleProducer/test/template_cfg.py \
    inputFile=my_g4.root \
    maxEvents=1000 \
    arm=1 \
    station=2 \
    rp=3 \
    verbosity=1
```

Output: `g4_reco_output.root`

## 3. Visualize Results

```bash
# From any directory with the Python libraries available:
python3 ~/PPS/visualize_pps_tracks.py \
    --g4-input ../pps2_simulation_prototype/test_output.root \
    --output-dir plots/
```

Produces: `plots/truth_kinematics.png`, `hit_distributions.png`, etc.

## 4. Validate Quickly (No CMSSW)

```bash
python3 ~/PPS/test_producer_simple.py
```

Checks:
- Event loading
- RecHit creation
- HepMC truth generation
- Hit efficiency

## Troubleshooting

### "No such file" for G4 ntuple
```bash
# Check file exists and has trees
root -b -e "auto f = TFile::Open('test_output.root'); f->ls();"
```

### Empty tracks in output
- Need ≥3 planes with hits per event
- Check hit positions: `-8 < x,y < 8` mm
- Run with `verbosity=2` to debug

### CMSSW build fails
```bash
# Clean and rebuild
scramv1 b clean SimPPS/G4NtupleProducer
scramv1 b SimPPS/G4NtupleProducer
```

### Python import errors
```bash
# Ensure uproot, numpy, matplotlib installed
pip install uproot numpy matplotlib
```

## Key Files

| File | Purpose |
|------|---------|
| `plugins/PPSG4NtupleToRecHitProducer.cc` | Main producer |
| `plugins/PPSG4NtupleValidator.cc` | Validation analyzer |
| `python/ppsG4NtupleToRecHitProducer_cfi.py` | Producer config |
| `python/matching_cff.py` | Reconstruction rewiring |
| `test/template_cfg.py` | Full reconstruction config |
| `test/make_test_ntuple.py` | Synthetic data generator |
| `README.md` | Full documentation |

## Configuration Parameters

**Producer settings** (`template_cfg.py` arguments):

```python
# Detector location
arm = 0              # 0=sector45, 1=sector56
station = 0          # 0=210m, 2=220m
rp = 3               # Roman Pot index

# Hit resolution (mm, sigma not variance)
sigmaHitX = 0.029    # 100um / sqrt(12)
sigmaHitY = 0.043    # 150um / sqrt(12)

# Control
verbosity = 0        # 0=quiet, 1=info, 2=debug
maxEvents = -1       # -1 = all events
```

## Output Collections

In `g4_reco_output.root`:

```
ctppsPixelLocalTracks       → Fitted tracks in each sensor
ctppsLocalTrackLiteProducer → Lite version (ready for proton fit)
ppsG4NtupleToRecHitProducer → RecHits + truth (HepMCProduct)
```

## One-Line Workflows

### Full pipeline (G4 → Reco → Plots)
```bash
cd ~/PPS/pps2_simulation_prototype && \
nix develop --command bash -c "./build/PPS2sim macros/run.mac" && \
cd ~/PPS/CMSSW_14_0_0 && \
cmsRun src/SimPPS/G4NtupleProducer/test/template_cfg.py inputFile=../pps2_simulation_prototype/test_output.root && \
python3 ~/PPS/visualize_pps_tracks.py --g4-input ../pps2_simulation_prototype/test_output.root --output-dir plots/
```

### Test with synthetic data
```bash
cd ~/PPS && \
python3 test_producer_simple.py && \
python3 visualize_pps_tracks.py --g4-input test_energy_angle_xy_10k.root --output-dir plots_syn/
```

## Performance Notes

- **~66% hit efficiency** (per plane, typical)
- **~0.1–0.5 ms/event** processing time
- **~100 MB/100k events** memory usage
- **10k events** from G4 → Reco → Plots: ~2 min total

## Need Help?

1. Check `README.md` (full documentation)
2. Run `test_producer_simple.py` (debug producer)
3. Increase `verbosity=2` in reconstruction config
4. Look at plot outputs for geometry issues

---

Last update: 30 July 2026
