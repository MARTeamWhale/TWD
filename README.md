# TWD
Toothed Whale Detector (tentative name)
Version 1.3

***NOTE: This project is under development!***

The Toothed Whale Detector (TWD) is a MATLAB-based software package designed to assist analysts in determining the presence of beaked and sperm whales from underwater acoustic recordings. Much of it is based on code that was developed by Simone Baumann-Pickering and Joy Stanistreet to analyze passive acoustic beaked whale data (Baumann-Pickering *et al.*, 2013, 2016; Stanistreet *et al.*, 2017). The sperm whale component also makes use of methods developed by Beslin *et al.* (2018).

Note that TWD itself is not a click detection and classification system. It is a tool to help analysts efficiently process large acoustic datasets to find periods of time where toothed whale clicks were present. While it does include some basic click discrimination parameters to facilitate this process, TWD does not attempt to produce accurate species-specific classification scores for clicks. The beaked whale component also does not include its own click detection routine.

TWD is currently divided into two modules:

- BWD (Beaked Whale Detector)
- SWD (Sperm Whale Detector)

Each module involves a two-step process:

1) Automatic detection of beaked or sperm whale acoustic events based on clicks
2) Interactive validation of detection events and (for beaked whales) determination of species present

The two steps are controlled via "master" script files located directly within the module directories (i.e., *_BWD* and *_SWD*). Use these scripts to set parameters and run the code.

More information is available under the *Documentation* folder (note that documentation specific to the latest version of TWD, version 1.3, is still being developed).


## References
- Baumann-Pickering, S. *et al.* (2013). Species-specific beaked whale echolocation signals. *J. Acoust. Soc. Am.* **134**, 2293–2301.

- Baumann-Pickering, S., Trickey, J. S., Wiggins, S. M. & Oleson, E. M. (2016). Odontocete occurrence in relation to changes in oceanography at a remote equatorial Pacific seamount. *Mar. Mammal Sci.* **32**, 805–825.

- Beslin, W. A. M., Whitehead, H. & Gero, S. (2018). Automatic acoustic estimation of sperm whale size distributions achieved through machine recognition of on-axis clicks. *J. Acoust. Soc. Am.* **144**, 3485–3495.

- Stanistreet, J. E. *et al.* (2017). Using passive acoustic monitoring to document the distribution of beaked whale species in the western North Atlantic Ocean. *Can. J. Fish. Aquat. Sci.* **74**, 2098–2109.
