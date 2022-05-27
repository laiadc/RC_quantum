# Adapting reservoir computing to solve the Schrödinger equation

Adaptation of Reservoir Computing to integrate the time-dependent Schrödinger equation, propagating an initial wavefunction in time. 

The results are illustrated in a set of [Jupyter](https://jupyter.org/) notebooks. The whole code is the result of the work in <a href = "https://arxiv.org/abs/" target="_blank"> this paper</a>. Any contribution or idea to continue the lines of the proposed work will be very welcome.

**Remark**: We recommend the readers to view the notebooks locally and in *Trusted* mode for nicer presentation and correct visualization of the figures. 

In this work, a Reservoir Computing-based model is develop to propagate quantum wavefunctions with time. Since such wavefunctions are complex-valued high-dimensional arrays the reservoir computing formalism needs to be extended to cope with complex-valued data.  Furthermore,  we propose a multi-step learning strategy that avoids overfitting the training data. 


<p align="center"><img src="https://github.com/laiadc/RC_quantum/blob/main/figures/MultiRC.PNG"  align=middle width=600pt />
</p>
<p align="center">
<em>Architecture of multi-step training of the reservoir computing model </em>
</p>

We illustrate the performance of our adapted reservoir computingmethod by application to four standard problems in molecular vibrational dynamics: the 1D harmonic oscillator, Morse Hamiltonian, polynomial potential and 2D harmonic oscillator.

## Notebooks

All the notebooks used for this work can be found inside the folder **notebooks** .

**Remark**: Some of the training data could not be uploaded because it exceeded the maximum size allowed by GitHub. The notebooks provide the code to obtain such training data. 

### [RC harmonic Oscillator 1D.ipynb](https://github.com/laiadc/RC_quantum/blob/main/notebooks/RC%20Harmonic%20Oscillator%201D.ipynb)
Application of the adapted Reservoir Computing model to solve the 1D harmonic oscillator.

### [RC Morse 1D.ipynb](https://github.com/laiadc/RC_quantum/blob/main/notebooks/RC%20Morse%201D.ipynb)
Application of the adapted Reservoir Computing model to solve the 1D Morse Hamiltonian.

### [RC polynomial 1D.ipynb](https://github.com/laiadc/RC_quantum/blob/main/notebooks/RC%20Random%20potential%201D.ipynb)
Application of the adapted Reservoir Computing model to solve the 1D polynomial potential.

### [RC harmonic Oscillator 2D.ipynb](https://github.com/laiadc/RC_quantum/blob/main/notebooks/RC%20Harmonic%20Oscillator%202D.ipynb)
Application of the adapted Reservoir Computing model to solve the 2D harmonic oscillator.

### [Figures.ipynb](https://github.com/laiadc/RC_quantum/blob/main/notebooks/Figures.ipynb)
This notebook summarizes the results of the paper and provides interactive plots created with *plotly* library.

### BibTex reference format for citation for the Code
```
@misc{RCQuantum,
title={Adapting reservoir computing to solve the Schrödinger equation},
url={https://github.com/laiadc/RC_quantum/},
note={Adaptation of Reservoir Computing to integrate the time-dependent Schrödinger equation, propagating an initial wavefunction in time.},
author={L. Domingo and J. Borondo and F. Borondo},
  year={2022}
}


