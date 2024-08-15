# Charged fluids near interfaces: Integral equation theory (Ph.D. Thesis)

## Web App

TL;DR

See [https://phd.databooth.com.au](https://phd.databooth.com.au) - a Streamlit web app deployed via Google Cloud Run.

## Background

I submitted my Ph.D. thesis `Charged fluids near interfaces: Integral equation theory` in 1998.
As part of the submission I included some electronic resources at the following link (no https in those days!):

[http://www.chem.usyd.edu.au/~boothm/thesis/thesis.html](http://www.chem.usyd.edu.au/~boothm/thesis/thesis.html)

Alas the University of Sydney (School of Chemistry) didn't seem to prioritise the longevity of these materials (I did ask others if they had a copy of my `FORTRAN` code, however to date nothing has materialised.) 

![Online Resources cover page](https://github.com/Mjboothaus/PhD-Thesis/blob/deb94173c8ee3b2591c19219e2ea813f33181ae1/docs/Thesis_OnlineResources_Cover.jpg)

## Bulk fluid correlation functions

The code for this was obtained from this project http://pyoz.vrbka.net which now seems to be discontinued.

and a fork is available here: https://github.com/ctk3b/pyoz.

I also found another (parallel) repo here: https://github.com/elvissoares/PyOZ (untested).

## TODOs

- This is still a work-in-progress - I think the LJ fluids work ok - but the charged fluids should have short-range correlations that approach zero (which is currently not the case).
- I have made an attempt at an "OO-refactoring" but I think this is only partial.
- The automatic wiring-in of the bulk fluid code would also be ideal.
- I think the bulk (and maybe singlet OZ) code can probably be simplified - I think it uses far from ideal discrete data structures.
- The capture of the solver text can possibly use a now in-built Streamlit widget but this is untested (can't find it as yet - maybe still to be released - 1.37 currently).
