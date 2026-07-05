# Blood Group Predictor

An interactive Shiny app that predicts ABO and Rh blood types independently, 
using genotype input, serology (antisera reaction) input, or parental cross probability.

🔗 **Live demo:** [View Here](https://m8wbzf-riffat-naz.shinyapps.io/blood_group_predictor/)

## Features
- ABO and Rh predicted independently — toggle Rh on/off without affecting ABO
- Three modes: Genotype, Serology, Parental Cross (Punnett square probabilities)
- Built with R Shiny

## Disclaimer
This tool is for educational purposes only. Rh factor is modeled as a single gene (D/d), whereas real-world Rh biology involves multiple genes and weak/partial variants. ABO logic assumes standard inheritance and doesn't account for rare exceptions (e.g. cis-AB, Bombay phenotype). Not intended for medical, diagnostic, paternity, or legal use.
