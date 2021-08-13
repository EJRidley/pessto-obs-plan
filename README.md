# ePESSTO+ Observation Planning

Useful code for planning ePESSTO+ observations.

**Dependencies: numpy, matplotlib, requests, astropy, astroplan**

**Requires ePESSTO+ Marshall login details in login.json**

### On First Run
core.py will generate:
* login.json (user should enter their details into this new file)
* ignore_list.txt
* graphs/
* outputs/

### Tips
ignore_list.txt can be used to remove targets which were planned for previous/current nights' observations, but have not yet been removed from the "classification targets" queue on the marshall.

The final lines of core.py produce the altitude plots. These can be adjusted to change observation date, choose observing location, filter the priorities of the events displayed, etc. The internal rank system corresponds to the priority system on the marshall. A helpful dictionary exists at the top of core.py which displays the relationships.

The angular distance to the moon is shown at the peak altitude of each target, and the moon's illumination for that night is displayed in each plot's title.
