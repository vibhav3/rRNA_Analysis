# rRNA Operon Neutral Evolution Analysis Tool

Current Version can be viewed at: https://r-rna-analysis.vercel.app/

## Overview

This tool tests whether rRNA operons in *E. coli* evolve neutrally or under selection by simulating neutral expectations and comparing them to experimental data from long-term evolution experiments.

**Three hypotheses:**
1. **Neutral evolution**: Mutation and gene conversion are balanced at equilibrium
2. **Purifying selection**: Selection maintains operon similarity
3. **Diversifying selection**: Selection favors operon differences

## How It Works

### Input Data
- CSV files with polymorphism data from sequencing timepoints
- Format: `Strain_Generation_FileID.csv` (e.g., `Ara-1_5000gen_762B.csv`)
- Contains: position, reference base, alternate base, coverage, allele frequency
- Compares all operons to rrnB (reference operon)

### Neutral Model Simulation

The tool simulates forward-time evolution incorporating:

**1. Mutation** (increases divergence)
```
Rate: 1×10⁻¹⁰ per bp per generation (normal strains)
      5×10⁻⁹ per bp per generation (mutator strains) (50 fold, lower end estimate)
Source: Long-Term Evolution Experiment (LTEE) estimates
```

**2. Gene Conversion** (decreases divergence)
```
Rate: 8.6×10⁻⁶ per operon per donor per generation
Tract length: ~50-250 bp
Source: Gifford et al. 2020
```

**3. Initial Divergence**
```
User-specified based on ancestral strain data
Typical value: ~0.002 (0.2% of sites differ)
```

### Algorithm

For each of 100+ simulation replicates:
```
1. Initialize 7 operons with starting differences
2. For each generation (0 to 50,000):
   - Apply mutations stochastically to each operon
   - Apply gene conversion events between operon pairs
   - Every 2,500 generations: measure divergence from rrnB
3. Calculate statistics: mean, median, 25th/75th percentiles
```

**Divergence metric**: Average fraction of sites differing from reference across all polymorphic positions

### Visualization

**Two synchronized plots:**
1. **Neutral expectation** (simulation results)
   - Mean trajectory (blue line)
   - Confidence bands: 25th-75th percentile (light blue)
   - Median (dashed line)

2. **Observed data** (experimental results)
   - Each strain plotted as separate trajectory
   - Mutator strains highlighted in red
   - Same Y-axis scale for direct comparison

## Interpreting Results

| Pattern | Interpretation | Biological Meaning |
|---------|---------------|-------------------|
| **Within neutral bands** | Consistent with neutral evolution | No evidence for selection |
| **Below neutral bands** | Less divergence than expected | Purifying selection or increased gene conversion |
| **Above neutral bands** | More divergence than expected | Diversifying selection or reduced gene conversion |
| **Flat trajectory at equilibrium** | Mutation-conversion balance | Neutral evolution with stable equilibrium |

## Key Model Insights

### Equilibrium Dynamics
Without gene conversion:
```
Divergence ∝ μt (linear increase)
```

With gene conversion:
```
Divergence → μL/γ (plateau at equilibrium)
where μ = mutation rate, L = length, γ = gene conversion rate
```

### Observed Data Analysis
```
For each timepoint and file, we calculate the frequency-weighted,
per-site dissimilarity from the rrnB reference by summing the frequencies of all
observed polymorphisms, grouped by operon. This sum is then normalized by
the total possible sites that could vary within each operon, rather than
the number of polymorphisms actually found. This approach produces a robust,
frequency-weighted measure of sequence divergence that is directly comparable to
neutral model expectations.
```

### Mutator Strain Predictions
- Higher equilibrium divergence (10-100× more mutations)
- Similar trajectory shape (gene conversion unchanged)
- Same time to equilibrium (if Ne large)

## Technology Stack

- **Frontend**: Next.js 14 (React 18)
- **Visualization**: Recharts
- **Styling**: Tailwind CSS
- **Language**: TypeScript
- **Icons**: Lucide React

## Installation
```bash
# Clone repository
git clone https://github.com/yourusername/rrna-analysis
cd rrna-analysis

# Install dependencies
npm install

# Run development server
npm run dev

# Open browser to http://localhost:3000
```

## Usage

1. **Upload CSV files**: Select multiple files from your experimental timepoints
2. **Set parameters**: 
   - Strain type (normal vs mutator)
   - Initial divergence (from ancestral data)
   - Number of simulations (100-1000)
3. **Run simulations**: Generate neutral expectations
4. **Compare**: Visual inspection + statistical testing

## File Structure
```
├── app/
│   ├── globals.css          # Global styles
│   ├── layout.tsx           # Root layout
│   └── page.tsx             # Home page
├── components/
│   └── RRNAAnalysis.tsx     # Main analysis component
├── package.json             # Dependencies
├── tailwind.config.js       # Tailwind configuration
└── tsconfig.json            # TypeScript configuration
```

## Statistical Testing (Future Development)

Current version provides visual comparison. Planned additions:

- **Goodness-of-fit tests**: Chi-square, Kolmogorov-Smirnov
- **Trajectory comparison**: Linear regression slope testing
- **Variance analysis**: Test for heterogeneity across strains
- **Bayesian parameter estimation**: ABC methods for rate inference

## Scientific Basis

### Key Assumptions
- Constant mutation rate (may vary with mutator phenotypes)
- Constant gene conversion rate (may change with recombination machinery)
- Large effective population size (drift negligible)
- No epistasis or linkage disequilibrium
- Independent evolution of each operon (except via gene conversion)

### Literature Support
- Mutation rates (ancestral): Wielgoss et al., 2011 
- Gene Converssion Rates: Gifford et al., 2020
- Concerted evolution: Liao 1999, Ganley & Kobayashi 2007

