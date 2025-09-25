# gge-analysis-macro
Usable in ImageJ to analyze particle size on PAGE scans: <br />
    Supports multiple standard lanes <br />
    Quartic regression

Also by default calculates amount of protein in different LDL bins (determined by MW): <br />
    Bins can be edited in macro code


Instructions for installation:
1. Download <a href='/particle-size-calculator/size_calculator.ijm' target='_blank'>the size calculator file</a> and return to this page by clicking the back button on your browser.
<img title='download instructions' alt='download instructions' src='/images/download_instructions.png'>
2. Open ImageJ.
3. In ImageJ, go to Plugins>Macros and select Install.
4. In the File Explorer Window, navigate to the location of size-calculator.ijm (probably in Downloads) and select it. Press Open.
5. The macro can now be initiated from Plugins>Macros>size-calculator


Instructions for use: <br>
IMPORTANT: Read all instructions carefully. If you receive an error, take note of the error output and restart from step 2.
1. Open ImageJ.
2. Open gel file via File>Open and navigate to gel image.
3. Initiate macro via Plugins>Macros>size-calculator. If the macro cannot be found, follow installation instructions above.
4. Enter gel identification number (Currently only used to title data output in Log).
5. Bound gel with rectangle, leaving a small border so as to not cut off any parts of the gel. You will have a chance later to more precisely mark the origin of each lane. Press OK.
6. For each standard lane:
    1. Enter number of standards.
    1. Enter standard weights.
    1. Position rectangle to capture entire length of lane. It can be quite thin. Press OK. Press 'Yes' if asked whether the lanes are horizontal.
    1. Read displayed instructions for marking origin and standard peaks and press OK.
    1. Select the origin. The origin will often present as a sharp peak signifying the edge of the gel. If you make a mistake, use Alt+Left Click to remove a point.
    1. Mark each standard peak and press the space bar once finished. Make sure the number of standard peaks you mark matches the number of standards you inputted for that lane.
7. For each sample lane:
    1. Position rectangle to capture entire length of lane. It can be quite thin. Press OK.
    1. Read displayed instructions for marking origin and unknowns and press OK.
    1. Select the origin.
    1. Select each unknown and press the space bar once finished.
    1. Read displayed instructions for marking baseline Y value and press OK.
    1. Select baseline Y value and press the space bar once finished.
    1. Continue analyzing lanes if desired.
8. The output in the log window will display the calculated MW/diameters of your selected peaks. By default, the macro will also try to calculate the amount of LDL in each LDL bin. Ignore this data if it is not relevant.

