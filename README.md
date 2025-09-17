# gge-analysis-macro
Usable in ImageJ to analyze particle size on PAGE scans:
    Supports multiple standard lanes
    Cubic linear regression

Also by default calculates amount of protein in different LDL bins (determined by MW):
    Bins can be edited in macro code


Instructions for use:
IMPORTANT: Read all instructions carefully.If you receive an error, take note of the error output and restart from step 2.
1. Open ImageJ.
2. Open gel file via File>Open and navigating to gel image.
3. Start macro via >Path> .Once path has been selected, macro can be run with Plugins>Tools>ImageJ GGE Macro
4. Enter gel X number (Currently only used to title data output in Log)
5. Bound gel with rectangle. You will have a chance later to more precisely mark the origin. Press OK.
6. For each standard lane:
    a. Enter number of standards
    b. Enter standard weights
    c. Position rectangle to capture entire length of lane. It can be quite thin. Press OK. Press 'Yes' if asked whether the lanes are horizontal.
    d. Read displayed instructions and press OK.
    e. Select the origin. The origin will often present as a sharp peak signifying the edge of the gel. If you make a mistake, use Alt+Left Click to remove a point.
    e. Mark each standard peak and press OK. Make sure the number of standard peaks you mark matches the number of standards you inputted for that lane.
7. For each sample lane:
    a. Position rectangle to capture entire length of lane. It can be quite thin.
    b. Mark origin of lane on plot and press OK. The origin will often present as a sharp peak signifying the edge of the gel. 
    c. Mark each unknown and press OK. 
8. The output in the log window will display the calculated MW/diameters of your selected peaks. By default, the macro will try to calculate the amount of LDL in each LDL bin. Ignore this data if it is not relevant.

