# Adapting reservoir computing to solve the Schrödinger equation

Adaptation of Reservoir Computing to integrate the time-dependent Schrödinger equation, propagating an initial wavefunction in time. 

The results are illustrated in a set of [Jupyter](https://jupyter.org/) notebooks. The whole code is the result of the work in <a href = "https://arxiv.org/abs/" target="_blank"> this paper</a>. Any contribution or idea to continue the lines of the proposed work will be very welcome.

**Remark**: We recommend the readers to view the notebooks locally and in *Trusted* mode for nicer presentation and correct visualization of the figures. 

In this work, a Reservoir Computing-based model is develop to propagate quantum wavefunctions with time. Since such wavefunctions are complex-valued high-dimensional arrays the reservoir computing formalism needs to be extended to cope with complex-valued data.  Furthermore,  we propose a multi-step learning strategy that avoids overfitting the training data. 


<p align="center"><img src="https://github.com/laiadc/RC_quantum/tree/main/figures/MultiRC.PNG"  align=middle width=600pt />
</p>
<p align="center">
<em>Architecture of multi-step training of the reservoir computing model </em>
</p>
