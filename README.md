# centauro

A python implementation of the [Centauro algorithm](https://arxiv.org/abs/2006.10751) (*by Miguel Arratia, Yiannis Makris, Duff Neill, Felix Ringer, and Nobuo Sato*) for jet clustering in the Breit-Frame of deep inelastic scattering (DIS). The algorithm consistently captures the jets initiate by the struck quark that is usually found in the very backward direction (*eta ~ -infinity*) where conventional algorithms such as  anti-KT (longitudinally invariant) fail. 

The package also includes implementation of the [event level grooming procedure](https://arxiv.org/abs/2101.02708) (*by Yiannis Makris*) using an extension of the modified-MassDrop Tagging (mMDT) procedure. The grooming is asymmetric favoring backward-energetic particles which intuitively result from the fragmentation of the struck quark in the leading-order DIS process.

For both implementation the algorithm assumes the convention that the target hadron moves along the +z direction. 

Example of input files with corresponding python scrips for importing the data are provided. 

A Mathematica notebook is also provided to visualization of the clustering and grooming procedure on the event. The visualization is performed on the unfolded sphere with final state particles are represented by discs with area representing their energy.

#### The main usage involves:
___
`jet[` *list_of_fourmomenta*, *daughter_1*, *daughter_2* `= 0` `]`

*list_of_fourmomenta* :  A python list of the jet constituents four-momenta of type *<float>* e.g., ` [    [E1, p1x, p1y, p1z],    [E2, p2x, p2y, p2z], ...  ]`  .

*daughter_1* :  integer > 0 , record number of daughter-1 of type *<int>*. For final state particles this has to b set to 0. 

*daughter_2* :  integer > 0 , record number of daughter-2 of type *<int>*.  For final state particles this has to b set to 0. 

**Initialize** an object `jet` based on the specs given.
___
`centauro_clustering = [` *list_of_jets*, *R* `=1` `]` 

*list_of_jets* : A list of object `jet` to be clustered according to the Centauro algorithm. 

*R* :  positive number for jet radius of type *<float>*.  Typical values are 0.1 - 1.5.

**Returns** a list of objects `jet`. Those are the inclusive jets according to Centauro algorithm.
___
`centauro_grooming = [` *list_of_jets*, *z_cut* `=0.1` `]` 

*list_of_jets* : A list of objects `jet` to be clustered based on the centauro measure and then groomed during the declustering process. Usually the input jets are the single-particle, final-state jets.

*z_cut* : number in the range [0, 1] as the grooming parameter.  Typical values are 0.01 - 0.5 .

**Returns** a list of objects `jet`. Those are the particles that pass grooming.
___

  
