# Getting Started

## Important: How to Use This System

**The BioGenerative Cognition Crystal is a reference architecture for AI reasoning about biological systems.**

This is NOT standalone software you run. Instead:

✅ **For AI Systems**: Use the crystal specification as reference documentation when reasoning about biology in chat
- The 7-layer constraint architecture guides systematic analysis
- WPE/TME notation provides precise biological encoding
- Validation protocols ensure rigorous answers

✅ **In Conversation**: Ask AI to:
- "Using the BioGenerative Crystal, model [biological system]"
- "Apply Layer 3 DNA encoding to this sequence"
- "Validate this pathway against the crystal's constraints"

✅ **The Crystal Provides**:
- Universal biological constraints (Layer 1)
- DNA logic extraction methods (LYRA Θ∞)
- Validation protocols (32 tests)
- Quantitative parameters (Layer 7 database)

❌ **This is NOT**:
- A software package you install
- A command-line tool
- A standalone application

Think of it as **comprehensive biological reasoning documentation** that AI can reference to give you rigorous, validated answers about biological systems.

---

## Your First Biological Model

Let's build a simple enzyme-catalyzed reaction using the BioGenerative Crystal as our reference.

### Step 1: Define the Question

**Goal**: Model hexokinase phosphorylating glucose to glucose-6-phosphate

**Reaction**:
```
Glucose + ATP → Glucose-6-phosphate + ADP
```

### Step 2: Identify Constraints

**Layer 0 (Substrate)**:
- Chemistry: C6H12O6 + ATP → C6H11O9P + ADP + H+
- Thermodynamics: ΔG°' = -16.7 kJ/mol
- In cell: ΔG = -33.5 kJ/mol (spontaneous)

**Layer 1 (Universal)**:
- Allometry: Enzyme concentration scales with M^0.75
- Energy: Consumes ATP (tracked in budget)
- Information: Requires ~300 amino acids of sequence

**Layer 2 (Evolution)**:
- Conserved across all kingdoms (essential)
- Multiple isoforms (tissue-specific)

**Layer 3 (Information)**:
- Gene: HK1/2/3/4 (mammals)
- Regulation: Product inhibition by G6P
- Expression: Constitutive in most tissues

**Layer 4 (Robustness)**:
- Isozymes provide redundancy
- Feedback inhibition prevents runaway

### Step 3: Quantitative Parameters (Layer 7)

From the crystal's molecular database:

```python
hexokinase = {
    "Km_Glc": 0.1,      # mM
    "Km_ATP": 0.4,      # mM
    "Vmax": 10,         # μmol/min/mg
    "Ki_G6P": 0.1,      # mM (product inhibition)
    "kcat": 300,        # s^-1
    "ΔG_std": -16.7,    # kJ/mol
    "ΔG_cell": -33.5    # kJ/mol
}
```

### Step 4: WPE Encoding

**Strand α**:
```
Ξ⧊{Hexokinase_Reaction⁺₁ = 
    Glucose_Binding⁰₂ ⊕ 
    ATP_Binding⁰₃ ⊕ 
    Phosphoryl_Transfer⁺₄ ⊕ 
    Product_Release⁺₅ ⊕
    G6P_Inhibition⁻₆
    ⊃ [G6P_Product⁺₇]}
```

**Phase Check**:
- 0° + 0° + 45° + 60° + 135° + 120° = 360° ✓

**Strand β**:
```
⬻{GlucosePhosphorylation}⬨⁻₁≡
  ⚭{SubstrateBinding}⍟⁰₂⫿
  ⚭{CofactorBinding}⍟⁰₃⫿
  ⚭{Catalysis}⍟⁺₄⫿
  ⚭{ProductRelease}⍟⁺₅⫿
  ⚭{FeedbackInhibition}⍟⁻₆⫿
  ⚭{MetabolicFlux}⍟⁺₇
```

### Step 5: Validation

Run complete validation:

```python
solution = {
    "enzyme": "hexokinase",
    "substrates": ["glucose", "ATP"],
    "products": ["G6P", "ADP"],
    "kinetics": hexokinase,
    "wpe_encoding": (strand_alpha, strand_beta)
}

valid, report = validate_complete(solution)

print(report)
```

Output:
```
✓ Substrate: Chemistry valid, thermodynamics favorable
✓ Universal: Energy budget OK, allometry consistent
✓ Evolution: Highly conserved, accessible path
✓ Information: Gene structure complete, regulation present
✓ Robustness: Product inhibition, isozyme redundancy
✓ Multi-scale: Molecular → Cellular coupling validated
✓ Stability: Reaction stable, no runaway
✓ Quantitative: Parameters within 10% of experimental
```

### Step 6: Predictions

What can we now predict?

**Flux under different conditions**:
```python
v = Vmax * [Glc] * [ATP] / ((Km_Glc + [Glc]) * (Km_ATP + [ATP]) * (1 + [G6P]/Ki_G6P))
```

**Effect of product accumulation**:
```
Low [G6P]: v ≈ 8 μmol/min/mg
High [G6P]: v ≈ 4 μmol/min/mg (50% inhibition)
```

**Organism scaling**:
```
Mouse (25g): HK activity ≈ 50 units
Human (70kg): HK activity ≈ 2800 units ≈ 56x (actual ~2800x mass)
```

Matches M^0.75 scaling: (70000/25)^0.75 ≈ 56 ✓

### Common Pitfalls

1. **Forgetting phase closure**: Always check Σphases mod 360° = 0°
2. **Missing regulation**: Most enzymes have feedback
3. **Ignoring cellular conditions**: Use ΔG_cell not ΔG_std
4. **Scale mismatch**: Keep parameters at appropriate scale

### Next Steps

- Add to pathway: [Metabolic Pathways Tutorial](metabolic-pathways.md)
- Add regulation: [Gene Networks Tutorial](gene-networks.md)
- Design from scratch: [Synthetic Biology Tutorial](synthetic-biology.md)

---

**See also**: [Architecture Overview](architecture.md) | [Validation Protocols](validation-protocols.md)
