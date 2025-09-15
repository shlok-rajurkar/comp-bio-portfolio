# gge-analysis-macro
Usable in ImageJ to analyze particle size on PAGE scans:
    Supports multiple standard lanes
    Cubic linear regression

Also by default calculates amount of protein in different LDL bins (determined by MW):
    Bins can be edited in macro code


Instructions for use:
IMPORTANT: When the macro prompts you to take certain actions, take those actions AND THEN press OK.
1. Open ImageJ.
2. Open gel file via File>Open and navigating to gel image.
3. Start macro via Plugins>Macros>Install>Krauss_Lab>Core A>GGE>Tools>_______ .Once path has been selected, macro can be run with Plugins>Tools>ImageJ GGE Macro
4. Enter gel X number (Currently only used to title data output in Log)
5. Bound gel with rectangle. You will have a chance later to more precisely mark the origin. Press OK.
