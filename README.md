# PMD-TMS
Fast And Accurate Computational E-field Dosimetry for Group-Level Transcranial Magnetic Stimulation Targeting.<br /><br />
Transcranial magnetic stimulation (TMS) is used to study brain function and treat mental health disorders. During TMS, a coil placed on the scalp induces an E-field in the brain that modulates its activity. TMS is known to stimulate regions that are exposed to a large E-field. Clinical TMS protocols prescribe a coil placement based on scalp landmarks. There are inter-individual variations in brain anatomy that result in variations in the TMS-induced E-field at the targeted region and its outcome. These variations across individuals could in principle be minimized by developing a large database of head subjects and determining scalp landmarks that maximize E-field at the targeted brain region while minimizing its variation using computational methods. However, this approach requires repeated execution of a computational method to determine the E-field induced in the brain for a large number of subjects and coil placements. We developed a probabilistic matrix decomposition-based approach for rapidly evaluating the E-field induced during TMS for a large number of coil placements due to a pre-defined coil model. Our approach can determine the E-field induced in over 1 Million coil placements in 9.5 hours, in contrast, to over 5 years using a brute-force approach. After the initial set-up stage, the E-field can be predicted over the whole brain within 2-3 milliseconds and to 2% accuracy. We tested our approach in over 200 subjects and achieved an error of < 2% in most and < 3.5% in all subjects. We will present several examples of bench-marking analysis for our tool in terms of accuracy and speed. Furthermore, we will show the methods’ applicability for group-level optimization of coil placement for illustration purposes only


## Authors
| Author | Affiliation | Email |
| --- | --- | --- |
| Nahian I. Hasan | Elmore Family School of Electrical and Computer Engineering, Purdue University, WL, USA | nahianhasan1994@gmail.com |
| Dezhi Wang | Elmore Family School of Electrical and Computer Engineering, Purdue University, WL, USA | wang5355@purdue.edu |
| Luis J. Gomez | Elmore Family School of Electrical and Computer Engineering, Purdue University, WL, USA | ljgomez@purdue.edu |


## Please cite the original paper as well as this repository as follows.
@article{hasan2023fast,<br />
  title={Fast And Accurate Population Level Transcranial Magnetic Stimulation via Low-Rank Probabilistic Matrix Decomposition (PMD)},<br />
  author={Hasan, Nahian Ibn and Wang, Dezhi and Gomez, Luis},<br />
  journal={bioRxiv},<br />
  pages={2023--02},<br />
  year={2023},<br />
  publisher={Cold Spring Harbor Laboratory}<br />
}
<br /><br /><br />
@software{Hasan_PMD-TMS_2023,<br />
author = {Hasan, Nahian I. and Gomez, Luis},<br />
month = oct,<br />
title = {{PMD-TMS}},<br />
url = {https://github.com/NahianHasan/PMD-TMS.git},<br />
version = {0.0.0.1},<br />
year = {2023}<br />
}
