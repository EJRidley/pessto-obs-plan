# ePESSTO+ Observation Planning

Useful code for planning ePESSTO+ observations.

**Dependencies: numpy, matplotlib, requests, astropy, astroplan**

**Requires ePESSTO+ Marshall login details in login.json**

### On First Run
core.py will generate:
* login.json
* ignore_list.txt
* graphs/
* outputs/

### Useful Tips
ignore_list.txt can be used to remove targets which were planned for previous/current nights' observations, but have not yet been removed from the "classification targets" queue on the marshall.

The final lines of core.py produce the altitude plots. Feel free to change/add/remove these lines to your liking. The internal rank system corresponds to the priority system on the marshall. A helpful dictionary exists at the top of core.py which displays the relationships.
