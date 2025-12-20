Understand the matlab codes boarbitrary and port them to julia, add a similar example in the documentaion

em3dHH Structure 

```txt
em3dHH/
├── boem3dHHroot/           # Root-finding solver (fsolve)
│   ├── boem3dHHroot.m      # Main driver script
│   ├── fDrHH.m             # Dispersion relation D(ω,k)
│   ├── funAn.m, funBn.m, funCn.m  # Perpendicular integrals
│   ├── funZl.m             # Plasma dispersion function Z_l
│   ├── funIn.m             # Normalization integrals I_n
│   └── bo.in               # Input parameters
│
└── boHH/                   # Matrix eigenvalue solver (all solutions)
    ├── bo_main.m           # Main entry point
    ├── modules/
    │   ├── bo_initialize.m     # Parameter setup (612 lines)
    │   ├── bo_em3d_matrix.m    # Build 3×3 dielectric tensor
    │   ├── bo_kernel.m         # Main loop, calls eig()
    │   ├── bo_output.m         # Save results
    │   ├── bo_plot_all.m       # Plot all dispersion surfaces
    │   └── bo_plot_select.m    # Plot selected branches
    ├── genfv/                  # Distribution function generation
    │   ├── gen_fv2d.m          # Generate f(v∥,v⊥) on grid
    │   ├── expand_fv2d.m       # Hermite expansion
    │   ├── funa0lm2alm.m       # Coefficient conversion
    │   └── hermiteH0.m         # Hermite polynomials
    ├── input/
    │   ├── bo.in               # Species parameters
    │   ├── bo_setup.m          # Solver configuration
    │   └── bo_knwn.m           # Normalization settings
    └── output/                 # Test cases
        ├── ringbeam_Umeda12/   # Ring beam instability
        ├── shell_Min15/        # Shell distribution
        ├── firehose_Astfalk17/ # Firehose instability
        └── ice_Irvine18/       # Ion cyclotron emission
```



```txt
┌─────────────────┐
│    bo.in        │  Species parameters
└────────┬────────┘
         ▼
┌─────────────────┐
│  bo_initialize  │  Compute ωp, ωc, vt, ρc, etc.
└────────┬────────┘
         ▼
┌─────────────────┐     ┌──────────────┐
│   genfv/        │────▶│  a_lm coefs  │  (if ifv=0)
│  expand_fv2d    │     └──────────────┘
└────────┬────────┘
         ▼
┌─────────────────┐
│   bo_kernel     │  Loop over k-space
└────────┬────────┘
         ▼
    ┌────┴────┐
    ▼         ▼
┌───────┐ ┌────────┐
│ eig() │ │fsolve()│  Matrix or root-finding
└───┬───┘ └───┬────┘
    └────┬────┘
         ▼
┌─────────────────┐
│   ω(k) curves   │  Dispersion relations
└─────────────────┘
```