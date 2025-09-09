## **Physicochemical property**

#### _<font color=red>Molecular Weight</font>_

Contain hydrogen atoms. Optimal:100~600, based on Drug-Like Soft rule.

#### _<font color=red>Volume</font>_

Van der Waals volume.

#### _<font color=red>Density</font>_

Density = MW / Volume

#### _<font color=red>nHA</font>_

Number of hydrogen bond acceptors. Sum of all O and N. Optimal: 0~12, based on Drug-Like Soft rule.

#### _<font color=red>nHD</font>_

Number of hydrogen bond donors. Sum of all OHs and NHs. Optimal:0~7, based on Drug-Like Soft rule.

#### _<font color=red>nRot</font>_

Number of rotatable bonds. In some situation Amide C-N bonds are not considered because of their high rotational energy barrier. Optimal:0~11, based on Drug-Like Soft rule.

#### _<font color=red>nRing</font>_

Number of rings. Smallest set of smallest rings. Optimal:0~6, based on Drug-Like Soft rule.

#### _<font color=red>MaxRing</font>_

Number of atoms in the biggest ring. Number of atoms involved in the biggest system ring. Optimal:0~18, based on Drug-Like Soft rule.

#### _<font color=red>nHet</font>_

Number of heteroatoms. Number of non-carbon atoms (hydrogens included). Optimal:1~15, based on Drug-Like Soft rule.

#### _<font color=red>fChar</font>_

Formal charge. Optimal:-4 ~4, based on Drug-Like Soft rule

#### _<font color=red>nRig</font>_

Number of rigid bonds. Number of non-flexible bonds, in opposite to rotatable bonds. Optimal:0~30, based on Drug-Like Soft rule.

#### _<font color=red>Flexibility</font>_

Flexibility = nRot / nRig

#### _<font color=red>Stereo Centers</font>_

Number of stereocenters. Optimal: ≤ 2, based on Lead-Like Soft rule.

#### _<font color=red>TPSA</font>_

Topological polar surface area. Sum of tabulated surface contributions of polar fragments. Optimal:0~140, based on Veber rule.

#### _<font color=red>logS</font>_

**Data Description**: The logarithm of aqueous solubility value. The first step in the drug absorption process is the disintegration of the tablet or capsule, followed by the dissolution of the active drug. Low solubility is detrimental to good and complete oral absorption, and early measurement of this property is of great importance in drug discovery.

**Results interpretation**: The predicted solubility of a compound is given as the logarithm of the molar concentration (log mol/L). Compounds in the range from -4 to 0.5 log mol/L will be considered proper.

**References**: admetlab2.0

#### _<font color=red>logP</font>_

**Data Description**: The logarithm of the n-octanol/water distribution coefficient. log P possess a leading position with considerable impact on both membrane permeability and hydrophobic binding to macromolecules, including the target receptor as well as other proteins like plasma proteins, transporters, or metabolizing enzymes.

**Results interpretation**: Results interpretation: The predicted logP of a compound is given as the logarithm of the molar concentration (log mol/L). Compounds in the range from 0 to 3 log mol/L will be considered proper.

**References**: admetlab2.0

#### _<font color=red>logD7.4</font>_ _(update)_

**Data Description**: The logarithm of the n-octanol/water distribution coefficients at pH=7.4. To exert a therapeutic effect, one drug must enter the blood circulation and then reach the site of action. Thus, an eligible drug usually needs to keep a balance between lipophilicity and hydrophilicity to dissolve in the body fluid and penetrate the biomembrane effectively. Therefore, it is important to estimate the n-octanol/water distribution coefficients at physiological pH (logD7.4) values for candidate compounds in the early stage of drug discovery.

**Results interpretation**: The predicted logD7.4 of a compound is given as the logarithm of the molar concentration (log mol/L). Compounds in the range from 1 to 3 log mol/L will be considered proper.

**References**: Duan YJ, Fu L, Zhang XC, et al. Improved GNNs for Log D7.4 Prediction by Transferring Knowledge from Low-Fidelity Data. J Chem Inf Model. 2023;63(8):2345-2359. doi:10.1021/acs.jcim.2c01564

#### _<font color=red>melting point</font>_ _(new)_

**Data Description**: Melting point (the temperature in °C at which a chemical in the solid state changes to a liquid state). The melting point for 9385 chemicals was obtained from the database in EPI Suite (USEPA 2009).

**Results interpretation**: The predicted melting point of a compound is expressed in degrees Celsius (°C). Melting points below 25°C are classified as liquids, while melting points above 25°C are classified as solids.

**References**: Dong J, Wang N N, Liu K Y, et al. Chemometrics and Intelligent Laboratory Systems, 2017, 171:65-73.

#### _<font color=red>normal boiling point</font>_ _(new)_

**Data Description**: The normal boiling point is defined as the temperature at which a chemical boils at atmospheric pressure. The data set for this endpoint was obtained from the boiling point data contained in EPI Suite (USEPA 2009). 41 chemicals were removed from the data set since they were previously shown to be badly predicted and had experimental values which were significantly different (>50K) from other sources such as NIST(NIST 2010) and LookChem (Lookchem.com 2011). The final data set contained 5759 chemicals.

**Results interpretation**: The predicted melting point of a compound is expressed in degrees Celsius (°C). A normal boiling point below 25°C is categorized as a gas.

**References**: Dong J, Wang N N, Liu K Y, et al. Chemometrics and Intelligent Laboratory Systems, 2017, 171:65-73.

#### _<font color=red>pka</font>_ _(new)_

**Data Description**: Acid-base dissociation constant (pKa) value represents the strength of a drug molecule's acidity or basicity. Understanding the pKa values of drug molecules is crucial for predicting their solubility, absorption, distribution, and metabolism in different physiological environments. The modeling data of pKa values are all macroscopic pKa values, and there is a mathematical relationship between micro-pKa and macro-pKa, as shown in equations:
$$
Acid:pK^{macro}_{a}=-Log( \sum^{N}_{i=1}10^{-pKa _{(i)} ^{micro}})
$$

$$
Base:pK^{macro}_{a}=Log( \sum^{N}_{i=1}10^{pKa _{(i)} ^{micro}})
$$
where Micros-pKa refers to the pKa value of specific dissociation sites, while macro-pKa represents the cumulative effect of multiple dissociation sites. When a compound contains only one dissociation site, microscopic and macroscopic pKa values are identical.

**Results interpretation**: **pka_acidic** denotes the pKa of acidic sites, with smaller values indicating stronger acidity. **pka_basic** signifies the pKa of basic sites, with larger values suggesting weaker acidity of the conjugate acid and stronger basicity of the corresponding site.

**References**: Wu, J., Wan, Y., Wu, Z., Zhang, S., Cao, D., Hsieh, C. Y., & Hou, T. (2022). MF-SuP-pKa: multi-fidelity modeling with subgraph pooling mechanism for pKa prediction. Acta Pharmaceutica Sinica B.

## **Medicinal Chemistry**

#### _<font color=red>QED</font>_

**Data Description**: A measure of drug-likeness based on the concept of desirability. QED is calculated by integrating the outputs of the desirability functions based on eight drug-likeness related properties, including MW, log P, NHBA, NHBD, PSA, Nrotb, the number of aromatic rings (NAr), and the number of alerts for undesirable functional groups. Here, average descriptor weights were used in the calculation of QED.

**Results interpretation**: The mean QED is 0.67 for the attractive compounds, 0.49 for the unattractive compounds and 0.34 for the unattractive compounds considered too complex.

**Empirical decision**: > 0.67: excellent (green); ≤ 0.67: poor (red)

**References**: Bickerton G R, Paolini G V, Besnard J, et al. Quantifying the chemical beauty of drugs[J]. Nat Chem, 2012, 4(2): 90-8.

#### _<font color=red>SAScore</font>_ _(update)_

**Data Description**: Synthetic accessibility score is designed to estimate ease of synthesis of drug-like molecules, based on a combination of fragment contributions and a complexity penalty. The score is between 1 (easy to make) and 10 (very difficult to make).

**Results interpretation**: ES: Easy to synthesize; HS: Hard to synthesize. The output value represents the probability of being easily synthesizable, ranging from 0 to 1.

**Empirical decision**: 1: excellent (green); 0: poor (red)

**References**: Ertl P, Schuffenhauer A. Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions[J]. J Cheminform, 2009, 1(1): 8.

#### _<font color=red>GASA</font>_ _(update)_

**Data Description**: Synthetic accessibility score is designed to estimate ease of synthesis of drug-like molecules, based on a combination of fragment contributions and a complexity penalty.

**Results interpretation**: ES: Easy to synthesize; HS: Hard to synthesize; The output value represents the probability of being difficult to synthesize, ranging from 0 to 1.

**Empirical decision**: 1: excellent (green); 0: poor (red)

**References**: Yu J, Wang J, Zhao H, et al. Organic Compound Synthetic Accessibility Prediction Based on the Graph Attention Mechanism. J Chem Inf Model. 2022;62(12):2973-2986. doi:10.1021/acs.jcim.2c00038

#### _<font color=red> Fsp3 </font>_

**Data Description**: Fsp3, the number of sp3 hybridized carbons/total carbon count, is used to determine the carbon saturation of molecules and characterize the complexity of the spatial structure of molecules. It has been demonstrated that the increased saturation measured by Fsp3 and the number of chiral centers in the molecule increase the clinical success rate, which might be related to the increased solubility, or the fact that the enhanced 3D features allow small molecules to occupy more target space.

**Results interpretation**: Fsp3 ≥ 0.42 is considered a suitable value.

**Empirical decision**: ≥ 0.42：excellent (green); ＜0.42: poor (red)

**References**: Lovering F, Bikker J, Humblet C. Escape from flatland: increasing saturation as an approach to improving clinical success[J]. J Med Chem, 2009, 52(21): 6752-6.

#### _<font color=red> MCE-18 </font>_

**Data Description**: MCE-18 stands for medicinal chemistry evolution in 2018, and this measure can effectively score molecules by novelty in terms of their cumulative sp3 complexity. It can effectively score structures by their novelty and current lead potential in contrast to simple and in many cases false positive sp3 index, and given by the following equation:

$$
MCE18 = (AR+NAR+CHRIRAL+SPIRO+\frac{sp^3+Cyc-Acyc}{1+sp^3}\times Q^1)
$$

where AR is the presence of an aromatic or heteroaromatic ring (0 or 1), NAR is the presence of an aliphatic or a heteroaliphatic ring (0 or 1), CHIRAL is the presence of a chiral center (0 or 1), SPIRO is the presence of a spiro point (0 or 1), sp3 is the portion of sp3-hybridized carbon atoms (from 0 to 1), Cyc is the portion of cyclic carbons that are sp3 hybridized (from 0 to 1), Acyc is a portion of acyclic carbon atoms that are sp3 hybridized (from 0 to 1), and Q1 is the normalized quadratic index.

**Results interpretation**: < 45: uninteresting, trivial, old scaffolds, low degree of 3D complexity and novelty; 45-63: sufficient novelty, basically follow the trends of currently observed in medicinal chemistry; 63-78: high structural similarity to the compounds disclosed in patent records; >78: need to be inspected visually to assess their target profile and drug-likeness.

**Empirical decision**: ≥ 45：excellent (green); ＜45: poor (red)

**References**: Ivanenkov Y A, Zagribelnyy B A, Aladinskiy V A. Are We Opening the Door to a New Era of Medicinal Chemistry or Being Collapsed to a Chemical Singularity?[J]. J Med Chem, 2019, 62(22): 10026-10043.

#### _<font color=red> NPscore </font>_

**Data Description**: The Natural Product-likeness score is a useful measure which can help to guide the design of new molecules toward interesting regions of chemical space which have been identified as “bioactive regions” by natural evolution. The calculation consists of molecule fragmentation, table lookup, and summation of fragment contributions.

**Results interpretation**: The calculated score is typically in the range from −5 to 5. The higher the score is, the higher the probability is that the molecule is a NP.

**References**: Ertl P, Roggo S, Schuffenhauer A. Natural product-likeness score and its application for prioritization of compound libraries[J]. J Chem Inf Model, 2008, 48(1): 68-74.

#### _<font color=red> Lipinski Rule </font>_

**Content**: MW≤500; logP≤5; Hacc≤10; Hdon≤5

**Results interpretation**: If two properties are out of range, a poor absorption or permeability is possible, one is acceptable.

**Empirical decision**:  < 2 violations：excellent (green)；≥2 violations: poor (red)

**References**: Lipinski C A, Lombardo F, Dominy B W, et al. Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings[J]. Adv Drug Deliv Rev, 2001, 46(1-3): 3-26.

#### _<font color=red>Pfizer Rule</font>_

**Content**: logP > 3; TPSA < 75

**Results interpretation**: Compounds with a high log P (>3) and low TPSA (<75) are likely to be toxic.

**Empirical decision**: two conditions satisfied: poor (red); otherwise: excellent (green)

**References**: Hughes J D, Blagg J, Price D A, et al. Physiochemical drug properties associated with in vivo toxicological outcomes[J]. Bioorg Med Chem Lett, 2008, 18(17): 4872-5.

#### _<font color=red>GSK Rule</font>_

**Content**: MW ≤ 400; logP ≤ 4

**Results interpretation**: Compounds satisfying the GSK rule may have a more favorable ADMET profile.

**Empirical decision**:  0 violations: excellent (green); otherwise: poor (red)

**References**:  Gleeson M P. Generation of a set of simple, interpretable ADMET rules of thumb[J]. J Med Chem, 2008, 51(4): 817-34.

#### _<font color=red>Golden Triangle </font>_

**Content**: 200 ≤MW ≤50; -2 ≤ logD ≤5

**Results interpretation**: Compounds satisfying the GoldenTriangle rule may have a more favourable ADMET profile.

**Empirical decision**: 0 violations: excellent (green); otherwise: poor (red)

**References**: Johnson T W, Dress K R, Edwards M. Using the Golden Triangle to optimize clearance and oral absorption[J]. Bioorg Med Chem Lett,2009,19(19):5560-4.

#### _<font color=red>PAINS</font>_

**Data Description**: Pan Assay Interference Compounds (PAINS) is one of the most famous frequent hitters filters, which comprises 480 substructures derived from the analysis of FHs determined by six target-based HTS assay. By application of these filters, it is easier to screen false positive hits and to flag suspicious compounds in screening databases. One of the most authoritative medicine magazines Journal of Medicinal Chemistry even requires authors to provide the screening results with the PAINS alerts of active compounds when submitting manuscripts.

**Results interpretation**: If the number of alerts is not zero, the users could check the substructures by the DETIAL button.

**References**: Baell J B, Holloway G A. New substructure filters for removal of pan assay interference compounds (PAINS) from screening libraries and for their exclusion in bioassays[J]. J Med Chem, 2010, 53(7): 2719-40.

#### _<font color=red>ALARM NMR Rule</font>_

**Data Description**: Thiol reactive compounds. There are 75 substructures in this endpoint.

**Results interpretation**: If the number of alerts is not zero, the users could check the substructures by the DETIAL button.

**References**: Huth J R, Mendoza R, Olejniczak E T, et al. ALARM NMR: a rapid and robust experimental method to detect reactive false positives in biochemical screens[J]. J Am Chem Soc, 2005, 127(1): 217-24.

#### _<font color=red>BMS Rule</font>_

**Data Description**: Undesirable, reactive compounds. There are 176 substructures in this endpoint.

**Results interpretation**: If the number of alerts is not zero, the users could check the substructures by the DETIAL button.

**References**: Pearce B C, Sofia M J, Good A C, et al. An empirical process for the design of high-throughput screening deck filters[J]. J Chem Inf Model, 2006, 46(3): 1060-8.

#### _<font color=red>Chelator Rule</font>_

**Data Description**: Chelating compounds. There are 55 substructures in this endpoint.

**Results interpretation**: If the number of alerts is not zero, the users could check the substructures by the DETIAL button.

**References**: Agrawal A, Johnson S L, Jacobsen J A, et al. Chelator fragment libraries for targeting metalloproteinases[J]. ChemMedChem, 2010, 5(2): 195-9.

#### _<font color=red>Colloidal aggregators</font>_ _(new)_

**Data Description**: Compound aggregation tends to start when the concentration is above the CAC and end as aggregators form with radius of approximately 30-600 nm. The resulting colloidal aggregators would non-specifically bind to the surface of proteins, thus inducing local protein unfolding, which usually results in destabilization or denaturation of enzymes

**Results interpretation**: Category 0: non-colloidal aggregators; Category 1: colloidal aggregators. The output value is the probability of being colloidal aggregators, within the range of 0 to 1.

#### _<font color=red>FLuc inhibitors</font>_ _(new)_

**Data Description**: Due to its unique catalysis mechanism, FLuc is widely used in a variety of HTS bioluminescence assays, especially in the assay which aims to study gene expression at the transcriptional level. However, the inhibition of Fluc by unexpected FLuc inhibitors would produce interference to HTS assays.

**Results interpretation**: Category 0: non-fLuc inhibitors; Category 1: fLuc inhibitors. The output value is the probability of being fLuc inhibitors, within the range of 0 to 1.

#### _<font color=red>Blue/Green fluorescence</font>_ _(new)_

**Data Description**: Fluorescence is the process by which a molecule, called fluorophore or fluorescent dye, absorbs a photon of light, exciting an electron to a higher energy state. Fluorophores have many applications, including as enzyme substrates, labels for biomolecules, cellular stains and environmental indicators. However, the appearance of fluorescent compound would produce interference to related HTS assays.

**Results interpretation**: Category 0: non-blue/green fluorescence; Category 1: blue/green fluorescence. The output value is the probability of being blue/green fluorescence, within the range of 0 to 1.

#### _<font color=red>Reactive compounds</font>_ _(new)_

**Data Description**: Chemical reactive compounds typically result in the chemical modification of reactive protein residues or, less frequently, the modification of nucleophilic assay reagents.

**Results interpretation**: Category 0: non-reactive compound; Category 1: reactive compound. The output value is the probability of being reactive compound, within the range of 0 to 1.

#### _<font color=red>Promiscuous compounds</font>_ _(new)_

**Data Description**: Promiscuous compounds refer to compounds that specifically bind to different macromolecular targets. These multiple interactions may include unintended targets, thus triggering adverse reactions and other safety issues.

**Results interpretation**: Category 0: non-promiscuous compound; Category 1: promiscuous compound. The output value is the probability of being promiscuous compound, within the range of 0 to 1.

## **Absorption**

#### _<font color=red>Caco-2 Permeability</font>_ _(update)_

**Data Description**: Before an oral drug reaches the systemic circulation, it must pass through intestinal cell membranes via passive diffusion, carrier-mediated uptake or active transport processes. The human colon adenocarcinoma cell lines (Caco-2), as an alternative approach for the human intestinal epithelium, has been commonly used to estimate in vivo drug permeability due to their morphological and functional similarities. Thus, Caco-2 cell permeability has also been an important index for an eligible candidate drug compound.

**Results interpretation**: The predicted Caco-2 permeability of a given compound is given as the log cm/s. A compound is considered to have a proper Cao-2 permeability if it has predicted value >-5.15log cm/s.

**Empirical decision**: > -5.15: excellent (green); otherwise: poor (red)

#### _<font color=red>PAMPA</font>_ _(new)_

**Data Description**: The Parallel Artificial Membrane Permeability Assay (PAMPA) is a method used to simulate the passive diffusion absorption of orally administered drugs in the gastrointestinal tract. This experiment employs the stirring double-sink PAMPA method to determine the permeability of compounds through PAMPA via passive diffusion. The resulting effective permeability (Peff) values are expressed in units of 10^-6 cm/s. By considering PAMPA in conjunction with Caco-2 and MDCK models, it is possible to overcome the limitations of individual models, thereby enabling a more accurate prediction of a drug's permeation and absorption characteristics in the body.

**Results interpretation**: The experimental data for Peff was logarithmically transformed (logPeff). Molecules with log Peff values below 2.0 were classified as low-permeability (Category 0), while those with log Peff values exceeding 2.5 were classified as high-permeability (Category 1).

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>MDCK Permeability</font>_

**Data Description**: Madin−Darby Canine Kidney cells (MDCK) have been developed as an in vitro model for permeability screening. Its apparent permeability coefficient, Papp, is widely considered to be the in vitro gold standard for assessing the uptake efficiency of chemicals into the body. Papp values of MDCK cell lines are also used to estimate the effect of the blood-brain barrier (BBB).

**Results interpretation**: The unit of predicted MDCK permeability is cm/s. A compound is considered to have a high passive MDCK permeability for a Papp > 20 x 10-6 cm/s, medium permeability for 2-20 x 10-6cm/s, low permeability for < 2 x 10-6cm/s.

**Empirical decision**:  >2 x 10-6cm/s: excellent (green), otherwise: poor (red)

#### _<font color=red>Pgp-inhibitor</font>_

**Data Description**: The inhibitor of P-glycoprotein. The P-glycoprotein, also known as MDR1 or 2 ABCB1, is a membrane protein member of the ATP-binding cassette (ABC) transporters superfamily. It is probably the most promiscuous efflux transporter, since it recognizes a number of structurally different and apparently unrelated xenobiotics; notably, many of them are also CYP3A4 substrates.

**Results interpretation**: Category 0: Non-inhibitor; Category 1: Inhibitor. The output value is the probability of being Pgp-inhibitor, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>Pgp-substrate</font>_

**Data Description**: As described in the Pgp-inhibitor section, modulation of P-glycoprotein mediated transport has significant pharmacokinetic implications for Pgp substrates, which may either be exploited for specific therapeutic advantages or result in contraindications.

**Results interpretation**: Category 0: Non-substrate; Category 1: substrate. The output value is the probability of being Pgp-substrate, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>HIA</font>_

**Data Description**: As described above, the human intestinal absorption of an oral drug is the essential prerequisite for its apparent efficacy. What’s more, the close relationship between oral bioavailability and intestinal absorption has also been proven and HIA can be seen an alternative indicator for oral bioavailability to some extent.

**Results interpretation**: A molecule with an absorbance of less than 30% is considered to be poorly absorbed. Accordingly, molecules with a HIA >30% were classified as HIA- (Category 0), while molecules with a HIA < 30% were classified as HIA+(Category 1). The output value is the probability of being HIA+, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>F20%</font>_

**Data Description**: For any drug administrated by the oral route, oral bioavailability is undoubtedly one of the most important pharmacokinetic parameters because it is the indicator of the eﬃciency of the drug delivery to the systemic circulation.

**Results interpretation**: Molecules with a bioavailability ≥ 20% were classified as F20%- (Category 0), while molecules with a bioavailability < 20% were classified as F20%+ (Category 1). The output value is the probability of being F20%+, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>F30%</font>_

**Data Description**:  For any drug administrated by the oral route, oral bioavailability is undoubtedly one of the most important pharmacokinetic parameters because it is the indicator of the eﬃciency of the drug delivery to the systemic circulation.

**Results interpretation**: Molecules with a bioavailability ≥ 30% were classified as F30%- (Category 0), while molecules with a bioavailability < 30% were classified as F30%+ (Category 1). The output value is the probability of being F30%+, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>F50%</font>_ _(new)_

**Data Description**:  For any drug administrated by the oral route, oral bioavailability is undoubtedly one of the most important pharmacokinetic parameters because it is the indicator of the eﬃciency of the drug delivery to the systemic circulation.

**Results interpretation**: Molecules with a bioavailability ≥ 50% were classified as F30%- (Category 0), while molecules with a bioavailability < 50% were classified as F50%+ (Category 1). The output value is the probability of being F50%+, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

## **Distribution**

#### _<font color=red>BCPR inhibitor</font>_ _(new)_

**Data Description**: Breast cancer resistance protein (BCRP), an ATP-binding cassette (ABC) efflux transporter, plays a critical role in multi-drug resistance (MDR) to anti-cancer drugs and drug–drug interactions. The prediction of BCRP inhibition can facilitate evaluating potential drug resistance and drug–drug interactions in early stage of drug discovery.

**Results interpretation**: Category 0: Non-inhibitor; Category 1: inhibitor. The output value is the probability of being inhibitor, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>OATP1B1 / 3 inhibitor</font>_ _(new)_

**Data Description**: The organic anion transporting polypeptide 1B1 (OATP1B1) and 1B3 (OATP1B3) are important hepatic uptake transporters. Inhibition of their normal functions could lead to drug-drug interactions.

**Results interpretation**: Category 0: Non-inhibitor; Category 1: inhibitor. The output value is the probability of being inhibitor, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>BSEP inhibitor</font>_ _(new)_

**Data Description**: Drug-induced cholestasis is a frequently observed side effect of drugs and is often caused by an unexpected interaction with the bile salt export pump (BSEP). BSEP is the key membrane transporter responsible for the transport of bile acids from hepatocytes into bile. Bile acid retention can be also caused by compounds that inhibit BSEP activity.

**Results interpretation**: Category 0: Non-inhibitor; Category 1: inhibitor. The output value is the probability of being inhibitor, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>MRP1 inhibitor</font>_ _(new)_

**Data Description**: Multidrug resistance protein 1 (MRP1), an integral transmembrane efflux transporter, belongs to the ATP-binding cassette (ABC) protein superfamily. MRP1 governs the absorption and disposition of a wide variety of endogenous and xenobiotic substrates including various drugs across organs and physiological barriers. Inhibition of their normal functions could lead to drug-drug interactions.

**Results interpretation**: Category 0: Non-inhibitor; Category 1: inhibitor. The output value is the probability of being inhibitor, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>PPB</font>_

**Data Description**: Plasma protein binding. One of the major mechanisms of drug uptake and distribution is through PPB, thus the binding of a drug to proteins in plasma has a strong influence on its pharmacodynamic behavior. PPB can directly influence the oral bioavailability because the free concentration of the drug is at stake when a drug binds to serum proteins in this process.

**Results interpretation**: A compound is considered to have a proper PPB if it has predicted value < 90%, and drugs with high protein-bound may have a low therapeutic index.

**Empirical decision**: ≤ 90%: excellent (green); otherwise: poor (red).

#### _<font color=red>VDss</font>_ _(update)_

**Data Description**:  The volume of distribution at steady state (VDss) is a fundamental pharmacokinetics (PK) property of drugs, which measures how effectively a drug molecule is distributed throughout the body. Along with the clearance (CL), it determines the half-life and, therefore, the drug dosing interval.

**Results interpretation**: The unit of predicted VDss is L/kg. A compound is considered to have a proper VDss if it has predicted VDss in the range of 0.04-20L/kg.

**Empirical decision**: 0.04-20: excellent (green); otherwise: poor (red)

#### _<font color=red>BBB Penetration</font>_

**Data Description**: Drugs that act in the CNS need to cross the blood–brain barrier (BBB) to reach their molecular target. By contrast, for drugs with a peripheral target, little or no BBB penetration might be required in order to avoid CNS side effects.

**Results interpretation**: The unit of BBB penetration is cm/s. Molecules with logBB > -1 were classified as BBB+ (Category 1), while molecules with logBB ≤ -1 were classified as BBB- (Category 0). The output value is the probability of being BBB+, within the range of 0 to 1.

**Empirical decision**:  0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>Fu</font>_

**Data Description**: The fraction unbound in plasms. Most drugs in plasma will exist in equilibrium between either an unbound state or bound to serum proteins. Efficacy of a given drug may be affect by the degree to which it binds proteins within blood, as the more that is bound the less efficiently it can traverse cellular membranes or diffuse.

**Results interpretation**: >20%: High Fu; 5-20%: medium Fu; <5% low Fu.

**Empirical decision**:  ≥ 5%: excellent (green);< 5%: poor (red).

## **Metabolism**

#### _<font color=red>CYP 1A2 / 2C19 / 2C9 / 2D6 / 3A4 / 2B6 (new) / 2C8 (new) inhibitor</font>_

#### _<font color=red>CYP 1A2 / 2C19 /2C9 / 2D6 / 3A4 / 2B6 (new) substrate</font>_

**Data Description**: Based on the chemical nature of biotransformation, the process of drug metabolism reactions can be divided into two broad categories: phase I (oxidative reactions) and phase II (conjugative reactions). There are 57 P450 isozymes in humans, and several of them serve as the major drug- metabolizing enzymes, including CYP1A2, CYP3A4, CYP2C9, CYP2C19, CYP2D6, CYP2B6, CYP2C8. Most of these CYPs responsible for phase I reactions are concentrated in the liver. In addition to the beneficial role in drug elimination, CYP-mediated metabolism could lead to the formation of toxic metabolites and occurrence of undesirable drug−drug interactions.

**Results interpretation**: Category 0: Non-substrate / Non-inhibitor; Category 1: substrate / inhibitor. The output value is the probability of being substrate / inhibitor, within the range of 0 to 1.

#### _<font color=red>HLM</font>_ _(new)_

**Data Description**: The human liver microsomal (HLM) stability assay is perhaps the one most commonly utilized for assessing clearance of chemicals by the human liver−the most important organ for drug metabolism.

**Results interpretation**: ≤ 30 min: unstable, > 30 min: stable. The unstable compounds were labeled as 1, and stable compounds were labeled as 0.

**Empirical decision**: 0-0.3:poor (red) ; 0.3-0.7: medium (yellow); 0.7-1.0(++): excellent (green)

## **Excretion**

#### <font color=red>$CL_{plasma}$</font>

**Data Description**: Plasma clearance ($CL_{plasma}$) is the most important pharmacokinetic parameter because it is the only one which controls the overall drug exposure (for a given bioavailability) and it is the parameter which allows computation of the dosage required to maintain an average steady-state plasma concentration.  $CL_{plasma}$ expresses the overall ability of the body to eliminate a drug by scaling the drug elimination rate (amount per time) by the corresponding plasma concentration level.

**Results interpretation**:  The unit of predicted  $CL_{plasma}$  penetration is ml/min/kg.  >15 ml/min/kg: high clearance; 5-15 ml/min/kg: moderate clearance; < 5 ml/min/kg: low clearance.

**Empirical decision**: 0-5: excellent (green); 5-15: medium (yellow); > 15 : poor (red)

#### <font color=red>$CL_{int}$</font> _(update)_

**Data Description**: Hepatic intrinsic clearance in the liver microsome ($CL_{int}$) is a parameter used to measure the metabolism and clearance rate of a drug within liver microsomes. $CL_{int}$ represents the metabolic and clearance rate specific to a particular drug within liver microsomes. It is an important pharmacokinetic parameter that helps to understand the drug's metabolism rate in the body, especially within the liver.

**Results interpretation**: The unit of predicted $CL_{int}$ is μL/min/mg. ≤ 5 μL/min/mg: low $CL_{int}$; 5-12 μL/min/mg: medium $CL_{int}$; ≥ 12 μL/min/mg: high $CL_{int}$.

**Empirical decision**: 0-5: excellent (green); 5-20: medium (yellow); > 20 :poor (red)

#### <font color=red>$T_{1/2}$</font> _(update)_

**Data Description**: The half-life ($T_{1/2}$) of a drug is a hybrid concept that involves clearance and volume of distribution, and it is arguably more appropriate to have reliable estimates of these two properties instead.

**Results interpretation**: The unit of predicted $T_{1/2}$ is n hours. ultra-short half-life drugs: $T_{1/2}$ < 1 hour;  short half-life drugs: $T_{1/2}$ between 1-4 hours; intermediate short half-life drugs: $T_{1/2}$ between 4-8 hours; long half-life drugs: $T_{1/2}$ > 8 hours.

## **Toxicology**

#### _<font color=red>DINeurot</font>_ _(new)_

**Data Description**: Drug-induced neurotoxicity (DINeurot) is a very common adverse reaction. Many clinical drugs and industrial chemicals can damage the nervous system in various ways. They can harm both the central nervous system and the peripheral nervous system, leading to neurological symptoms and episodes resembling psychosis. The extent of these damages to the nervous system may be temporary and reversible, or it may result in long-term, irreversible organic injuries.

**Results interpretation**: Category 0: non-neurotoxic (-); Category 1: neurotoxic (+). The output value is the probability of being neurotoxic (+), within the range of 0 to 1.

**Empirical decision**:  0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>Ototoxicity</font>_ _(new)_

**Data Description**: Ototoxicants are compounds that have the potential to harm the inner ear by either damaging the ear's structures directly or the nervous system.

**Results interpretation**: Category 0: non-ototoxicity (-); Category 1:  ototoxicity (+). The output value is the probability of being ototoxicity (+), within the range of 0 to 1.

**Empirical decision**:  0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>Hematotoxicity</font>_ _(new)_

**Data Description**: Chemical-induced hematotoxicity refers to adverse effects of chemicals on blood-forming organs such as bone marrow or the constituents of blood, including platelets, leuko- cytes and erythrocytes.

**Results interpretation**: Category 0: non-hematotoxicity (-); Category 1:  hematotoxicity (+). The output value is the probability of being hematotoxicity (+), within the range of 0 to 1.

**Empirical decision**:  0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>DIN</font>_ _(new)_

**Data Description**: Nephrotoxicity refers to the harmful effects that occur in the kidneys due to chemicals and medicines, known as nephrotoxicants, often resulting in their rapid deterioration. The kidneys are uniquely susceptible to drug-induced injury due to their high cardiac output and their role in the excretion of waste compounds from the body. Due to their pivotal role in concentrating and reabsorbing the glomerular filtrate, the kidney proximal tubular cells are particularly prone to elevated levels of circulating toxicants. Drug-induced nephrotoxicity (DIN) has been identified as a major contributor to both acute kidney injury and chronic kidney disease.

**Results interpretation**: Category 0: non-nephrotoxic (-); Category 1:  nephrotoxic (+). The output value is the probability of being nephrotoxic (+), within the range of 0 to 1.

**Empirical decision**:  0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>Genotoxicity</font>_ _(new)_

**Data Description**: Genotoxicity refers to the ability of harmful substances to damage genetic information in cells. Being exposed to chemical and biological agents can result in genomic instabilities and/or epigenetic alterations, which translate into a variety of diseases, cancer included.

**Results interpretation**: Category 0: non-Genotoxicity (-); Category 1:  Genotoxicity (+). The output value is the probability of being ototoxicity (+), within the range of 0 to 1.

**Empirical decision**:  0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>RPMI-8226</font>_ _(new)_

**Data Description**: The evaluation of compound cytotoxicity is a crucial step in the drug discovery process. The RPMI-8226 cell line is a type of multiple myeloma cell line, and understanding the effects of compounds on RPMI-8226 cells helps determine the toxicity of these compounds to this cell line.

**Results interpretation**: Category 0: non-cytotoxicity (-); Category 1:  cytotoxicity (+). The output value is the probability of being ototoxicity (+), within the range of 0 to 1.

**Empirical decision**:  0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>A549</font>_ _(new)_

**Data Description**: The evaluation of compound cytotoxicity is a critical step in the drug discovery process. A549 is a human non-small cell lung cancer cell line, and understanding the effects of compounds on A549 cells helps determine the toxicity of these compounds to this cell line.

**Results interpretation**: Category 0: non-cytotoxicity (-); Category 1:  cytotoxicity (+). The output value is the probability of being ototoxicity (+), within the range of 0 to 1.

**Empirical decision**:  0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>Hek293</font>_ _(new)_

**Data Description**: HEK293 cells are human embryonic kidney cells, and these cells are relatively sensitive to the external microenvironment. This sensitivity can be utilized for drug discovery and cell cytotoxicity assay.

**Results interpretation**: Category 0: non-cytotoxicity (-); Category 1:  cytotoxicity (+). The output value is the probability of being ototoxicity (+), within the range of 0 to 1.

**Empirical decision**:  0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>hERG inhibitor_cls10</font>_ _(new)_

**Data Description**: The human ether-a-go-go related gene. The During cardiac depolarization and repolarization, a voltage-gated potassium channel encoded by hERG plays a major role in the regulation of the exchange of cardiac action potential and resting potential. The hERG blockade may cause long QT syndrome (LQTS), arrhythmia, and Torsade de Pointes (TdP), which lead to palpitations, fainting, or even sudden death.

**Results interpretation**: Molecules with IC50 ≤10 μM are classified as hERG+ (Category 1), and molecules with IC50 > 10μM are classified as hERG- (Category 0). The output value is the probability of being hERG+, within the range of 0 to 1.

**Empirical decision**:  0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>hERG inhibitor</font>_

**Data Description**: The human ether-a-go-go related gene. The During cardiac depolarization and repolarization, a voltage-gated potassium channel encoded by hERG plays a major role in the regulation of the exchange of cardiac action potential and resting potential. The hERG blockade may cause long QT syndrome (LQTS), arrhythmia, and Torsade de Pointes (TdP), which lead to palpitations, fainting, or even sudden death.

**Results interpretation**: Molecules with IC50 ≤10μM or ≥50% inhibition at 10 μM were classified as hERG+ (Category 1), while molecules with IC50 >10μM or <50%  inhibition at 10μM were classified as hERG - (Category 0). The output value is the probability of being hERG+, within the range of 0 to 1.

**Empirical decision**:  0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>DINeurot</font>_ _(new)_

**Data Description**: Drug-induced neurotoxicity (DINeurot) is a very common adverse reaction. Many clinical drugs and industrial chemicals can damage the nervous system in various ways. They can harm both the central nervous system and the peripheral nervous system, leading to neurological symptoms and episodes resembling psychosis. The extent of these damages to the nervous system may be temporary and reversible, or it may result in long-term, irreversible organic injuries.

**Results interpretation**: Category 0: non-neurotoxic (-); Category 1: neurotoxic (+). The output value is the probability of being neurotoxic (+), within the range of 0 to 1.

**Empirical decision**:  0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)




#### _<font color=red>H-HT</font>_

**Data Description**: The human hepatotoxicity. Drug induced liver injury is of great concern for patient safety and a major cause for drug withdrawal from the market. Adverse hepatic effects in clinical trials often lead to a late and costly termination of drug development programs.

**Results interpretation**: Category 0: H-HT negative(-); Category 1: H-HT positive(+). The output value is the probability of being toxic, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)


#### _<font color=red>Hepatotoxicity</font>_

**Data Description**: Drug induced liver injury is of great concern for patient safety and a major cause for drug withdrawal from the market. Adverse hepatic effects in clinical trials often lead to a late and costly termination of drug development programs. This dataset was collected from three publicly available databases (comprehensive PharmaPendium, ChemIDplus, and manual toxicological literature searches).

**source**: https://pubs.acs.org/doi/10.1021/acs.chemrestox.2c00411

**Results interpretation**: Category 0: non-hepatotoxic; Category 1: hepatotoxic.
Note: These three datasets on hepatotoxicity—liver injury, Hepatotoxicity I, and Hepatotoxicity II—differ in their data sources.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>DILI</font>_

**Data Description**:  Drug-induced liver injury (DILI) has become the most common safety problem of drug withdrawal from the market over the past 50 years.Compounds collected from Zhang et al., (18) Ai et al., (22) Liew et al., (34) and Kotsampasakou et al. (36) were merged and duplicates were removed to build our development set. Compounds collected from Chen et al. (35) (DILIrank benchmark dataset) were used as our independent test set. For Chen et al.’s dataset, only compounds termed “most-DILI” and termed “non-DILI” were selected as positive and negative samples, respectively. Compounds in the independent test set that appeared in the development set were removed to ensure that there was no overlapping between the two datasets. To check for structural validity, all compounds were cross-referenced to PubChem to discard unidentified compounds. 

**Data source**: https://pubs.acs.org/doi/10.1021/acsomega.0c03866

**Results interpretation**:  Category 0: DILI negative(-); Category 1: DILI positive(+). The output value is the probability of being toxic, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>AMES Toxicity</font>_

**Data Description**: The Ames test for mutagenicity. The mutagenic effect has a close relationship with the carcinogenicity, and it is the most widely used assay for testing the mutagenicity of compounds.
**Results interpretation**: Category 0: AMES negative(-); Category 1: AMES positive(+). The output value is the probability of being toxic, within the range of 0 to 1.
**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>Rat Oral Acute Toxicity</font>_

**Data Description**: Determination of acute toxicity in mammals (e.g. rats or mice) is one of the most important tasks for the safety evaluation of drug candidates.

**Results interpretation**: Category 0: low-toxicity, > 500 mg/kg; Category 1: high-toxicity; < 500 mg/kg. The output value is the probability of being toxic, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>FDAMDD</font>_

**Data Description**: The maximum recommended daily dose provides an estimate of the toxic dose threshold of chemicals in humans.

**Results interpretation**: Category 1: FDAMDD positive(+), ≤ 0.011 mmol/kg -bw/day; Category 0: FDAMDD negative(-), > 0.011 mmol/kg-bw/day. The output value is the probability of being toxic, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>Skin Sensitization</font>_

**Data Description**: Skin sensitization is a potential adverse effect for dermally applied products. The evaluation of whether a compound, that may encounter the skin, can induce allergic contact dermatitis is an important safety concern.

**Results interpretation**: Category 1: Sensitizer; Category 0: Non-sensitizer. The output value is the probability of being toxic, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>Carcinogencity</font>_ _(update)_

**Data Description**: Among various toxicological endpoints of chemical substances, carcinogenicity is of great concern because of its serious effects on human health. The carcinogenic mechanism of chemicals may be due to their ability to damage the genome or disrupt cellular metabolic processes. Many approved drugs have been identified as carcinogens in humans or animals and have been withdrawn from the market.

**Results interpretation**: Category 1: carcinogens; Category 0: non-carcinogens. Chemicals are labelled as active (carcinogens) or inactive (non-carcinogens) according to their TD50 values. The output value is the probability of being toxic, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>Eye Corrosion / Irritation</font>_

**Data Description**: Assessing the eye irritation/corrosion (EI/EC) potential of a chemical is a necessary component of risk assessment. Cornea and conjunctiva tissues comprise the anterior surface of the eye, and hence cornea and conjunctiva tissues are directly exposed to the air and easily suffer injury by chemicals. There are several substances, such as chemicals used in manufacturing, agriculture and warfare, ocular pharmaceuticals, cosmetic products, and household products, that can cause EI or EC.

**Results interpretation**: Category 1: corrosives / irritants chemicals; Category 0: non-corrosives / non-irritants chemicals. The output value is the probability of being toxic, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>Respiratory Toxicity</font>_

**Data Description**: Among these safety issues, respiratory toxicity has become the main cause of drug withdrawal. Drug-induced respiratory toxicity is usually underdiagnosed because it may not have distinct early signs or symptoms in common medications and can occur with significant morbidity and mortality.Therefore, careful surveillance and treatment of respiratory toxicity is of great importance.

**Results interpretation**: Category 1: respiratory toxicants; Category 0: non-respiratory toxicants. The output value is the probability of being toxic, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>Bioconcentration Factor</font>_

**Data Description**: The bioconcentration factor BCF is defined as the ratio of the chemical concentration in biota as a result of absorption via the respiratory surface to that in water at steady state. It is used for considering secondary poisoning potential and assessing risks to human health via the food chain. The unit of BCF is log10(L/kg).

#### _<font color=red>IGC50</font>_

**Data Description**: 48 hour Tetrahymena pyriformis IGC50 (concentration of the test chemical in water in mg/L that causes 50% growth inhibition to Tetrahymena pyriformis after 48 hours). The unit of IGC50 is −log10[(mg/L)/(1000*MW)].

#### _<font color=red>LC50FM</font>_

**Data Description**: 96 hour fathead minnow LC50 (concentration of the test chemical in water in mg/L that causes 50% of fathead minnow to die after 96 hours). The unit of LC50FM is −log10[(mg/L)/(1000*MW)].

#### _<font color=red>LC50DM</font>_

**Data Description**: 48 hour Daphnia magna LC50 (concentration of the test chemical in water in mg/L that causes 50% of Daphnia magna to die after 48 hours). The unit of LC50DM is −log10[(mg/L)/(1000*MW)].

#### _<font color=red>NR-AR</font>_

**Data Description**: Androgen receptor (AR), a nuclear hormone receptor, plays a critical role in AR-dependent prostate cancer and other androgen related diseases. Endocrine disrupting chemicals (EDCs) and their interactions with steroid hormone receptors like AR may cause disruption of normal endocrine function as well as interfere with metabolic homeostasis, reproduction, developmental and behavioral functions.

**Results interpretation**: Category 1: actives ; Category 0: inactives. The output value is the probability of being AR agonists, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>NR-AR-LBD</font>_

**Data Description**: Androgen receptor (AR), a nuclear hormone receptor, plays a critical role in AR-dependent prostate cancer and other androgen related diseases. Endocrine disrupting chemicals (EDCs) and their interactions with steroid hormone receptors like AR may cause disruption of normal endocrine function as well as interfere with metabolic homeostasis, reproduction, developmental and behavioral functions.

**Results interpretation**: Category 1: actives ; Category 0: inactives. Molecules that labeled 1 in this bioassay may bind to the LBD of androgen receptor. The output value is the probability of being actives, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>NR-AhR</font>_

**Data Description**: The Aryl hydrocarbon Receptor (AhR), a member of the family of basic helix-loop-helix transcription factors, is crucial to adaptive responses to environmental changes. AhR mediates cellular responses to environmental pollutants such as aromatic hydrocarbons through induction of phase I and II enzymes but also interacts with other nuclear receptor signaling pathways.

**Results interpretation**: Category 1: actives ; Category 0: inactives. Molecules that labeled 1 may activate the aryl hydrocarbon receptor signaling pathway. The output value is the probability of being actives, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>NR-Aromatase</font>_

**Data Description**: Endocrine disrupting chemicals (EDCs) interfere with the biosynthesis and normal functions of steroid hormones including estrogen and androgen in the body. Aromatase catalyzes the conversion of androgen to estrogen and plays a key role in maintaining the androgen and estrogen balance in many of the EDC-sensitive organs.

**Results interpretation**: Category 1: actives ; Category 0: inactives. Molecules that labeled 1 are regarded as aromatase inhibitors that could affect the balance between androgen and estrogen. The output value is the probability of being actives, within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>NR-ER</font>_

**Data Description**: Estrogen receptor (ER), a nuclear hormone receptor, plays an important role in development, metabolic homeostasis and reproduction. Endocrine disrupting chemicals (EDCs) and their interactions with steroid hormone receptors like ER causes disruption of normal endocrine function. Therefore, it is important to understand the effect of environmental chemicals on the ER signaling pathway.

**Results interpretation**: Category 1: actives ; Category 0: inactives. The output value is the probability of being actives within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>NR-ER-LBD</font>_

**Data Description**: Estrogen receptor (ER), a nuclear hormone receptor, plays an important role in development, metabolic homeostasis and reproduction. Two subtypes of ER, ER-alpha and ER-beta have similar expression patterns with some uniqueness in both types. Endocrine disrupting chemicals (EDCs) and their interactions with steroid hormone receptors like ER causes disruption of normal endocrine function.

**Results interpretation**:  Category 1: actives ; Category 0: inactives. The output value is the probability of being actives within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>NR-PPAR-gamma</font>_

**Data Description**: The peroxisome proliferator-activated receptors (PPARs) are lipid-activated transcription factors of the nuclear receptor superfamily with three distinct subtypes namely PPAR alpha, PPAR delta (also called PPAR beta) and PPAR gamma (PPARg). All these subtypes heterodimerize with Retinoid X receptor (RXR) and these heterodimers regulate transcription of various genes. PPAR-gamma receptor (glitazone receptor) is involved in the regulation of glucose and lipid metabolism.

**Results interpretation**: Category 1: actives ; Category 0: inactives. The output value is the probability of being actives within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>SR-ARE</font>_

**Data Description**: Oxidative stress has been implicated in the pathogenesis of a variety of diseases ranging from cancer to neurodegeneration. The antioxidant response element (ARE) signaling pathway plays an important role in the amelioration of oxidative stress. The CellSensor ARE-bla HepG2 cell line (Invitrogen) can be used for analyzing the Nrf2/antioxidant response signaling pathway. Nrf2 (NF-E2-related factor 2) and Nrf1 are transcription factors that bind to AREs and activate these genes.

**Results interpretation**:  Category 1: actives ; Category 0: inactives. The output value is the probability of being actives within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>SR-ATAD5</font>_

**Data Description**: ATPase family AAA domain-containing protein 5. As cancer cells divide rapidly and during every cell division they need to duplicate their genome by DNA replication. The failure to do so results in the cancer cell death. Based on this concept, many chemotherapeutic agents were developed but have limitations such as low efficacy and severe side effects etc. Enhanced Level of Genome Instability Gene 1 (ELG1; human ATAD5) protein levels increase in response to various types of DNA damage.

**Results interpretation**: Category 1: actives ; Category 0: inactives. The output value is the probability of being actives within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>SR-HSE</font>_

**Data Description**: Heat shock factor response element. Various chemicals, environmental and physiological stress conditions may lead to the activation of heat shock response/ unfolded protein response (HSR/UPR). There are three heat shock transcription factors (HSFs) (HSF-1, -2, and -4) mediating transcriptional regulation of the human HSR.

**Results interpretation**: Category 1: actives ; Category 0: inactives. The output value is the probability of being actives within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>SR-MMP</font>_

**Data Description**: Mitochondrial membrane potential (MMP), one of the parameters for mitochondrial function, is generated by mitochondrial electron transport chain that creates an electrochemical gradient by a series of redox reactions. This gradient drives the synthesis of ATP, a crucial molecule for various cellular processes. Measuring MMP in living cells is commonly used to assess the effect of chemicals on mitochondrial function; decreases in MMP can be detected using lipophilic cationic fluorescent dyes.

**Results interpretation**: Category 1: actives ; Category 0: inactives. The output value is the probability of being actives within the range of 0 to 1.

**Empirical decision**: 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>SR-p53</font>_

**Data Description**: p53, a tumor suppressor protein, is activated following cellular insult, including DNA damage and other cellular stresses. The activation of p53 regulates cell fate by inducing DNA repair, cell cycle arrest, apoptosis, or cellular senescence. The activation of p53, therefore, is a good indicator of DNA damage and other cellular stresses.

**Results interpretation**: Category 1: actives ; Category 0: inactives. The output value is the probability of being actives within the range of 0 to 1.

**Empirical decision**:  0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0(++): poor (red)

#### _<font color=red>Acute Toxicity Rule</font>_

**Data Description**: Molecules containing these substructures may cause acute toxicity during oral administration. There are 20 substructures in this endpoint.

**Results interpretation**: If the number of alerts is not zero, the users could check the substructures by the DETIAL button.

#### _<font color=red>Genotoxic Carcinogenicity Rule</font>_

**Data Description**: Molecules containing these substructures may cause carcinogenicity or mutagenicity through genotoxic mechanisms.There are 117 substructures in this endpoint.

**Results interpretation**: If the number of alerts is not zero, the users could check the substructures by the DETIAL button.

#### _<font color=red>NonGenotoxic Carcinogenicity Rule</font>_

**Data Description**: Molecules containing these substructures may cause carcinogenicity through nongenotoxic mechanisms. There are 23 substructures in this endpoint.

**Results interpretation**: If the number of alerts is not zero, the users could check the substructures by the DETIAL button.

#### _<font color=red>Skin Sensitization Rule</font>_

**Data Description**: Molecules containing these substructures may cause skin irritation.There are 155 substructures in this endpoint. Molecules containing these substructures may cause skin irritation.

**Results interpretation**: If the number of alerts is not zero, the users could check the substructures by the DETIAL button.

#### _<font color=red>Aquatic Toxicity Rule</font>_

**Data Description**: Molecules containing these substructures may cause toxicity to liquid(water). There are 99 substructures in this endpoint.

**Results interpretation**: If the number of alerts is not zero, the users could check the substructures by the DETIAL button.

#### _<font color=red>NonBiodegradable Rule</font>_

**Data Description**: Molecules containing these substructures may be non-biodegradable. There are 19 substructures in this endpoint.

**Results interpretation**: If the number of alerts is not zero, the users could check the substructures by the DETIAL button.

#### _<font color=red>SureChEMBL Rule</font>_

**Data Description**: Molecules matching one or more structural alerts are considered to have MedChem unfriendly status. There are 164 substructures in this endpoint.
**Results interpretation**: If the number of alerts is not zero, the users could check the substructures by the DETIAL button.
