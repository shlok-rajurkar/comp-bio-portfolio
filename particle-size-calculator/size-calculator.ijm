//GGE particle size macro for ImageJ to replace GGE7A3 on NIHImage.

main();


function main() {
    initialize();
    getStandards();
    cubicRegression(stdRfVals, stdWeights);
    getLanes();
    }


function getStandards() {
    stdWeights = newArray(0);
    stdRfVals = newArray(0);
    moreStandards = true;
    while (moreStandards) {
        stdWeightsTempAndStdRfValsTemp = setStandards();
        arrayLength = stdWeightsTempAndStdRfValsTemp.length;
        stdWeightsToAdd = Array.slice(stdWeightsTempAndStdRfValsTemp, 0, (arrayLength/2));
        stdRfValsToAdd = Array.slice(stdWeightsTempAndStdRfValsTemp, (arrayLength/2), arrayLength);
        //Array.print(stdWeightsToAdd);
        //Array.print(stdRfValsToAdd);
        stdWeights = Array.concat(stdWeights, stdWeightsToAdd);
        stdRfVals = Array.concat(stdRfVals, stdRfValsToAdd);
        stdWeights = Array.sort(stdWeights);
        stdWeights = Array.reverse(stdWeights);
        stdRfVals = Array.sort(stdRfVals);
        //Array.print(stdWeights);
        //Array.print(stdRfVals);
        moreStandards = getBoolean("Add more standard lanes?");
    }
}

function getLanes() {
    moreLanes = true;
    while (moreLanes) {
        quantifyLane();
        moreLanes = getBoolean("Analyze more lanes?");
    }
}
var stdWeights;
var stdRfVals;
var slopeAndInterceptArray;
var imageWidth;
var imageHeight;
var xnum;
var c = 0;
var z = 0;
var t = 0;
var cubicCoeffArray;
var croppedGelWindow;
var originXVal;
var laneLength;

function calcRfVals(xVals, imgLength) {
    result = newArray(xVals.length);
    for (i = 0; i < xVals.length; i ++) {
        result[i] = xVals[i]/imgLength;
    }
    return result;
}

function log10Array(array) {
    logValues = newArray(array.length);
    for (i = 0; i < array.length; i++){
        logValues[i] = log(array[i])/log(10);
    }
    return logValues;
}

function inverseLog10Array(array) {
    inverseLogValues = newArray(array.length);
    for (i = 0; i < array.length; i++){
        inverseLogValues[i] = pow(10, array[i]);
    }
    return inverseLogValues;
}

function getRfValsFromLane() {
    run("Plot Lanes");

    run("Remove Overlay");

    setTool("multi-point");

    waitForUser("Mark origin of lane, then press OK.");

    getSelectionCoordinates(originXVal, originYVal);

    waitForUser("Mark the peaks with the multipoint tool, then press OK. \nOnly mark them from left to right, in increasing weight/diameter.");

    getSelectionCoordinates(xpoints, ypoints);

    nPoints = xpoints.length;

    xVals = newArray(nPoints);

    for (i = 0; i < nPoints; i++) {
        xVals[i] = xpoints[i];
    }

    xVals = Array.slice(xVals, 1, xVals.length);

    getDimensions(plotWidth, plotHeight, c, z, t);

    RfVals = calcRfVals(xVals, plotWidth-originXVal[0]);

    return RfVals;    
}

function OLSRegression(xValsForRegression, yValsForRegression) {
    xLength = xValsForRegression.length;
    yLength = yValsForRegression.length;

    xTotal = 0;
    yTotal = 0;

    for(i = 0; i < xLength; i++) {
        xTotal += xValsForRegression[i];
        yTotal += yValsForRegression[i];
    }
    xMean = xTotal/xLength;
    yMean = yTotal/yLength;

    numerator = 0;
    denominator = 0;

    for(i = 0; i < xLength; i++) {
        dx = xValsForRegression[i] - xMean;
        dy = yValsForRegression[i] - yMean;
        numerator += dx * dy;
        denominator += dx * dx;
    }

    regressionSlope = numerator/denominator;

    regressionIntercept = yMean - regressionSlope * xMean;

    return newArray(regressionSlope, regressionIntercept);
}

function calcStdCurveLinear(stdRfValues, log10StdWeights) {
    slopeAndInterceptArray = OLSRegression(stdRfValues, log10StdWeights);  
    return slopeAndInterceptArray;  
}

function calcStdCurveCubic(xVals, yVals) {
    n = xVals.length;

    sumX = newArray(7);
    for (i = 0; i <= 6; i++) sumX[i] = 0;
    sumXY = newArray(4);
    for (i = 0; i <= 3; i++) sumXY[i] = 0;

    for (i = 0; i < n; i++) {
        x = xVals[i];
        y = yVals[i];
        powX = newArray(7);
        powX[0] = 1;
        for (j = 1; j <= 6; j++) powX[j] = powX[j - 1] * x;

        for (j = 0; j <= 6; j++) sumX[j] += powX[j];
        for (j = 0; j <= 3; j++) sumXY[j] += powX[j] * y;
    }

    // Flattened 4x4 matrix A
    A = newArray(16);
    for (i = 0; i <= 3; i++) {
        for (j = 0; j <= 3; j++) {
            A[i * 4 + j] = sumX[i + j];
        }
    }

    B = sumXY;

    coeffs = gaussJordanFlat4x4(A, B);
    return coeffs; // [a, b, c, d]
}

function gaussJordanFlat4x4(A, B) {
    n = 4;

    for (i = 0; i < n; i++) {
        // Find non-zero pivot
        if (A[i * 4 + i] == 0) {
            for (j = i + 1; j < n; j++) {
                if (A[j * 4 + i] != 0) {
                    // Swap rows in A
                    for (k = 0; k < n; k++) {
                        temp = A[i * 4 + k];
                        A[i * 4 + k] = A[j * 4 + k];
                        A[j * 4 + k] = temp;
                    }
                    // Swap B
                    tmpB = B[i];
                    B[i] = B[j];
                    B[j] = tmpB;
                    break;
                }
            }
        }

        // Normalize row
        factor = A[i * 4 + i];
        for (j = 0; j < n; j++) A[i * 4 + j] /= factor;
        B[i] /= factor;

        // Eliminate other rows
        for (j = 0; j < n; j++) {
            if (j != i) {
                factor = A[j * 4 + i];
                for (k = 0; k < n; k++) {
                    A[j * 4 + k] -= factor * A[i * 4 + k];
                }
                B[j] -= factor * B[i];
            }
        }
    }

    return B;
}



function initialize() {


    //stdWeights = newArray(0);
    //stdRfVals = newArray(0);
    //slopeAndInterceptArray = newArray(0);
    imageWidth = 0
    imageHeight = 0
    xnum = ""

    run("Set Measurements...", "invert redirect=None decimal=1");
    run("Gel Analyzer Options...", "vertical=1 horizontal=1 label");

    xnum = getString("enter X number", "");
    
    setTool("Rectangle");

    getDimensions(imageWidth, imageHeight, c, z, t);
    makeRectangle(0.3*imageWidth, 0.18*imageHeight, 0.36*imageWidth, 0.53*imageHeight);
    waitForUser("Adjust rectangle to bound lanes and leave minimal space between top of image and top of gel. \nPressing OK will crop gel and rotate 90 degrees left.");

    run("Crop");
    wait(100);
    run("Rotate 90 Degrees Left");
    croppedGelWindow = File.name;
    print(xnum + " Analysis:");
}


function setStandards() {
    selectWindow(croppedGelWindow);
    numberOfStandards = getNumber("enter number of standards", 5);
    stdWeightsTemp = newArray(numberOfStandards);
    stdRfValsTemp = newArray(numberOfStandards);
    for(i = 0; i < numberOfStandards; i++){
        stdIndexDisplay = i + 1;
        stdWeightsTemp[i] = getNumber("enter weight of standard " + stdIndexDisplay, 0);
    }
    print("Std weights:");
    //Array.print(stdWeightsTemp);

    setTool("Rectangle");

    getDimensions(standardLaneWidth, standardLaneHeight, c, z, t);

    makeRectangle(0, (standardLaneHeight/2)-10 , standardLaneWidth, 20);
    waitForUser("Adjust rectangle to span lane with standards. \nIt can be quite thin as long as it contains some part of the lane.");

    //usePresetRectangle = getBoolean('Use ROI preset?');

    //if (usePresetRectangle) {
    //    makeRectangle(0, (imageHeight/2)-10 , imageWidth, 20);
    //    waitForUser("Adjust rectangle to span lane with standards. \nIt can be quite thin as long as it contains some part of the lane.");
    //}
    //else {
    //    waitForUser("Draw rectangle from left to right of lane with standards. \nIt can be quite thin as long as it contains some part of the lane. Pressing OK will select this lane.");
    //}

    
    run("Select First Lane");

    stdRfValsTemp = getRfValsFromLane();
    //Array.print(stdRfValsTemp);
    stdWeightsTempAndStdRfValsTemp = Array.concat(stdWeightsTemp, stdRfValsTemp);
    return stdWeightsTempAndStdRfValsTemp;

    
}

function cubicRegression(stdRfValues, stdWeights) {
    log10StdWeights = log10Array(stdWeights);
    cubicCoeffArray = calcStdCurveCubic(stdRfVals, log10StdWeights);
    //Array.print(cubicCoeffArray);
    return cubicCoeffArray;
}

function quantifyLane() {
    if (stdRfVals.length == 0 || stdWeights.length == 0){
        exit("Standards not set.");
    }

    selectWindow(croppedGelWindow);

    waitForUser("Drag rectangle to sample lane and press OK.");

    run("Select First Lane");

    laneRfVals = getRfValsFromLane();

    laneMolecularWeightsCalc = newArray(laneRfVals.length);
    for (i = 0; i < laneRfVals.length; i++) {
        x = laneRfVals[i];
        laneMolecularWeightsCalc[i] = cubicCoeffArray[0] + cubicCoeffArray[1]*x + cubicCoeffArray[2]*x*x + cubicCoeffArray[3]*x*x*x;
    }
    displayValueArray = inverseLog10Array(laneMolecularWeightsCalc);
    print("Calculated Diameters:");
    Array.print(displayValueArray);

    quantBins();

}

function inverseLog10Array(array) {
    inverseLogValues = newArray(array.length);
    for (i = 0; i < array.length; i++){
        inverseLogValues[i] = pow(10, array[i]);
    }
    return inverseLogValues;
}
//Test cases for XG4513/14 where stds are 100, 80, 60, 40, 20
//imageWidth = 1506;
//originXVal = 0;
//cubicCoeffArray = newArray(2.0490, 0.0341, -0.8621, -0.1735);
print("-------------------------------------------------------");
function calcPxFromBins() {
    bins = newArray(
    375, 339, 321, 315, 309, 303, 297, 291, 285, 272, 265, 256, 247, 242, 233, 220
    );
    bins = newArray(100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5);
    getDimensions(samplePlotWidth, samplePlotHeight, c, z, t);
    laneLength = samplePlotWidth-originXVal[0];
    binCount = bins.length; 
    everyPixel = newArray(laneLength);
    for (i = 0; i < laneLength; i ++) {
        everyPixel[i] = i;
    }
    //print('every pixel:');
    //Array.print(everyPixel);
    everyRfValue = newArray(laneLength);
    everyLogMW = newArray(laneLength);
    for (i = 0; i < laneLength; i ++) {
        everyRfValue[i] = everyPixel[i]/laneLength;
    }
    //print('every Rf value:');
    //Array.print(everyRfValue);
    for (i = 0; i < laneLength; i ++) {
            everyLogMW[i] = cubicCoeffArray[0] + cubicCoeffArray[1]*everyRfValue[i] + cubicCoeffArray[2]*everyRfValue[i]*everyRfValue[i] + cubicCoeffArray[3]*everyRfValue[i]*everyRfValue[i]*everyRfValue[i];
    }
    //print('everyMW');
   
    everyMW = inverseLog10Array(everyLogMW);
    //Array.print(everyMW);
    binPxValues = newArray(binCount);
    //Array.print(binPxValues);
    for (j = 0; j < bins.length; j ++) {
        //print("original MWs");
        //Array.print(everyMW);
        everyMWCopy = Array.copy(everyMW);
        //everyMWCopy = everyMW;
        // print("before subtraction:");
        // Array.print(everyMWCopy);
        for (i = 0; i < laneLength; i ++) {
            everyMWCopy[i] = abs(everyMWCopy[i] - bins[j]);
            
        }
        // print("current bin:");
        // print(bins[j]);
        // print("after subtraction:");
        // Array.print(everyMWCopy);

        Array.getStatistics(everyMWCopy, min, max, mean, stdDev);
        // print("min=" + min);

        for (h = 0; h < laneLength; h ++) {
            if (abs(everyMWCopy[h] - min)<1e-6) {
                binPxValues[j] = h;
                //print(binPxValues[j]);
                //print(j);
            }
        }
    }
    //Array.print(binPxValues);
    return binPxValues;
}



function quantBins() {
    binPxValues = calcPxFromBins();
    //binPxValues = newArray(380, 450, 512, 569, 623, 676, 727, 778, 828, 879, 932, 986, 1042, 1102, 1167, 1238, 1319, 1415, 1505, 1505);
    baselineY = getBackgroundConc();
    //selectWindow("Plots of " + croppedGelWindow)
    //yVals = newArray(549.3548, 608.3548, 429.7742, 304.8710, 255.2258, 228.2581, 220.6452, 221.4194, 218.3871, 219.3871, 220.7742, 219.3548, 223.3226, 216.1935, 213.6774, 213.3226, 208.7097, 201.9032, 198.6452, 190.3871, 184.8710, 179.1935, 173.2581, 169.7742, 174.5806, 179.2581, 181, 181.5806, 172.6129, 166.7419, 160.2581, 149.1935, 153.3226, 150.3226, 141.2258, 140.4839, 142.7097, 144.1935, 144.0323, 139.4516, 137.2581, 137.1935, 136.4194, 137.7419, 138.0968, 138.9355, 137.0645, 138.8065, 136.1935, 135.9032, 136.4194, 136.8710, 135.7742, 136.1613, 133.7097, 133.2903, 133.1935, 135.0968, 133.2581, 132.9677, 132.7097, 133, 134.6452, 133.7742, 134.1290, 134.1290, 135.3548, 134.2903, 134.0645, 134.2258, 134.7419, 134.5806, 135.9355, 136.6774, 137.1290, 137.3548, 136.9355, 136.7419, 137.2903, 137.1613, 137.3548, 137.5161, 138.9032, 139.1613, 138.8065, 140.2258, 144.2258, 150.8710, 141.5161, 140.9032, 141.7742, 140.8387, 141.1613, 142.0323, 142.2903, 143, 143.2581, 143.6774, 143.7419, 144.1290, 144.9032, 145.9677, 145.3871, 145.1613, 145, 145.1613, 145.9355, 145.9677, 146.5161, 146.1935, 146.0968, 146.5484, 146.9355, 146.5484, 146.7419, 147.5484, 148.3226, 148.3226, 147.5806, 148.3548, 150.8710, 149.0645, 148.1290, 149.1290, 148.9677, 148.8710, 148.9677, 149.0968, 148.6129, 149.5806, 150.9677, 150.5484, 150.3548, 150.1613, 149.8710, 151.0323, 150.2903, 150.0968, 151.7097, 151.4516, 151.2581, 152.4839, 152.5484, 152.3226, 151.8710, 152.7097, 154.8387, 153.2581, 153.4516, 153.1290, 153.3226, 153.3548, 153.7419, 154.2258, 156.0645, 157.9355, 161.0323, 155.6129, 155.9032, 156.4839, 158.4516, 157.6774, 157.4839, 158.1935, 158.3871, 159.2258, 161.4839, 164.7742, 160.4839, 160.4839, 160.0645, 160.3226, 161.4839, 162.7419, 160.7742, 161.5161, 161, 161.5484, 160.3226, 160.7742, 161.4839, 165.4516, 163.0323, 162.0323, 162.6129, 162.7097, 162.9355, 162.8065, 163.0323, 162.7097, 163.0323, 163.1290, 163.5484, 163.6452, 164.1290, 165.3548, 164.6452, 167.0968, 164.1935, 162.8710, 164.2903, 163.8387, 164.1613, 164.7097, 165.5484, 166.5484, 167.9677, 167.1935, 166.8065, 167.5484, 166.6774, 166.2258, 166.9355, 166.8387, 167.5806, 167.8387, 168.0968, 167.5161, 168.6129, 168.5484, 169.2258, 168.8387, 168.4516, 168.7419, 168.0968, 169.6452, 169, 169.8710, 168.4194, 167.9677, 168.1935, 169.2258, 168.2258, 168.5161, 169, 168.2581, 167.8710, 169.1290, 170.8387, 171.3871, 173.8065, 174.3871, 172.3871, 171.8387, 172.3548, 172.9677, 173.3548, 173.9355, 174.7097, 175.9677, 175.8387, 175.3871, 174.4516, 174.0968, 174.3871, 174.1290, 174.3226, 174.7097, 174.2258, 173.7419, 173.3548, 174.1290, 172.3871, 171.6452, 170.2581, 169.1290, 169.6452, 169.2581, 169.5161, 169.8387, 170.0645, 171, 170.9677, 171.7097, 170.7419, 171.1613, 173.1613, 173.2258, 172.2581, 172.6452, 173.1935, 173.7097, 174.2903, 178.1290, 174.7419, 175.7419, 177.3548, 178.7742, 180.9355, 180.5806, 180.2903, 179.2258, 180.3226, 180.6774, 180.9032, 181.4839, 182.2258, 183.2903, 185.1935, 185.3548, 185.1290, 185.6452, 186.6452, 188.5484, 190.0968, 191.6774, 193.1935, 194.9032, 195.3548, 196.3226, 198.7097, 199.8387, 200.4194, 202.0968, 203.4194, 204.6129, 206.2258, 207.7419, 209.1935, 211.6774, 213.5161, 215.5161, 218.0323, 221, 223.1290, 225.8387, 228.0968, 230.4194, 235.3548, 236.6129, 239.5484, 244.8387, 247, 250.8387, 255.7097, 259.3226, 266.3548, 269.5806, 273.5161, 280.0323, 286.1935, 292.1613, 296.4516, 304.3871, 310.0323, 317.8387, 326.0645, 332.1613, 340.8710, 350.3871, 358.0645, 370.0968, 382.5806, 393.3548, 404.7097, 414, 424.3548, 430.3226, 435.8065, 439.9032, 443.9032, 448.9355, 452.6452, 456.4839, 462.4839, 466.5806, 469.5806, 472.9677, 474.1290, 473.7742, 471.5484, 466.5806, 460.7419, 454.3226, 443.1935, 429.8065, 413.1613, 401.0645, 384.2581, 363.7419, 348.4194, 330.9032, 315.3548, 303.1935, 288.0968, 276, 267.7097, 257.0645, 248.2258, 243, 235.2258, 229.0645, 227.7742, 221.3548, 216.9032, 214.5161, 212.1935, 209.7097, 208.4516, 206.4839, 204.7742, 204.2903, 203.7419, 203.4516, 203.4194, 202.3548, 202.4194, 201.2903, 200.5806, 200.0645, 200.0323, 203.2581, 201, 199.6129, 199.8710, 199.0645, 199.1613, 198.5484, 198.1613, 197.9032, 200.7097, 202.0968, 200.3226, 199.0645, 200.0323, 199.2903, 199.2581, 198.6774, 198.8710, 198.8065, 199.4516, 198.9677, 199.4516, 200.4194, 200.3871, 201.0645, 200.8065, 200.9677, 200.7419, 200.7097, 200.3548, 199.4516, 199.7419, 199.6774, 200.4194, 199.3871, 199.9355, 199.7742, 199.4194, 199.5806, 199.6129, 199.3226, 199.5484, 199.3871, 199.1613, 199.5484, 198.8710, 201.3871, 200.2258, 198.8387, 199.0645, 199.0968, 198.4516, 198.3548, 198.2258, 198.1613, 198.2903, 198.4194, 198.0645, 198.1613, 197.9677, 198.1935, 197.9355, 198.0645, 199.3871, 199, 198.7419, 199.6452, 199.3226, 199.8387, 200.2581, 200.1935, 200.8065, 199.9677, 200.3871, 200.6774, 201.3871, 200.9355, 200.7742, 201.6452, 201.9032, 201.8387, 202.6129, 203.5806, 202.6774, 202.8065, 203.6774, 203.8387, 204.7742, 205.1613, 205.5161, 206.7097, 207.3226, 208.0645, 208.6774, 209.6774, 212.5806, 216.2581, 211.6129, 212.5806, 213.6774, 214.3871, 214.3548, 216.7097, 217.6452, 217.9355, 218.8387, 218.2581, 220, 221.9355, 221.3548, 222.4194, 223.0968, 223.9355, 224.3548, 226.8710, 228.2581, 228.8710, 232.4516, 232.4516, 234.1613, 235.7097, 238.5161, 238.5484, 243.6452, 263.0968, 265.1935, 260.6129, 254.4194, 251.3548, 253.8387, 256.5806, 258.2903, 260.8065, 264.7419, 266.0968, 270.1290, 274.2581, 277.0645, 281.0645, 285.4839, 288.8065, 293.7097, 298.4839, 301.9355, 306.9032, 311.8710, 316.6774, 321.7097, 329.7419, 335.0645, 342.0645, 349.9032, 358.4194, 371.6129, 390.0968, 409.0323, 442.4516, 491.1613, 536.9677, 604.4839, 692.6452, 766.9355, 862.5806, 958.6452, 1025.3226, 1099.8387, 1158.0323, 1187.4839, 1202, 1194.3226, 1170.7742, 1123.0968, 1053.4194, 992.0968, 915.9032, 837.1290, 777.3226, 703.2903, 629.9355, 581.4516, 518.6129, 458.0968, 418.0323, 371.5484, 328.8065, 303.1613, 273.6774, 250.4516, 235.4516, 219.9355, 207.3871, 200.9032, 193.1290, 185.6452, 181.1290, 177.3548, 173.0645, 170.0323, 167.0968, 164.4839, 162.2581, 160.5806, 158.3226, 157.2581, 156.3548, 154.8710, 154.1613, 153.8710, 152.6129, 152.4194, 152.4839, 152.5161, 152.7419, 151.4516, 151.4839, 151.7419, 152.3548, 151.2258, 150.8710, 150.3871, 150.5161, 150.8710, 150.4194, 150.9677, 152.0323, 150.1290, 149.3871, 149.2581, 149.0968, 148.7419, 148.9355, 149.5484, 148.8065, 148.4839, 148.7742, 149.1613, 149.3871, 148.9355, 148.8387, 148.8710, 150.4516, 153.4194, 150.9677, 150.3548, 150.7742, 149.5806, 149.5806, 149.2581, 149.1290, 149.9677, 149.9355, 151.3548, 153.6129, 150.8710, 151.6452, 150.6452, 148.5161, 148.6129, 148.8387, 148.4194, 148.4516, 148.9032, 149.2581, 151.1290, 151.4839, 149.4516, 148.4516, 148.6774, 149.6129, 149.6452, 149.3226, 148.4839, 148.6129, 148.7097, 148.8387, 148.5806, 148.4839, 152.1613, 148.7097, 150.8710, 147.5484, 149, 148.6774, 147.3226, 147.7097, 147.9032, 147.0968, 146.3226, 145.0323, 146.4839, 147.2258, 148.7097, 148.0323, 148.1613, 146.4839, 146.0323, 145.5806, 145.2903, 146.2258, 145.7419, 145.5806, 145.8387, 145.5806, 147.3871, 148.2903, 146.0968, 146.9355, 148.7419, 147.0323, 147.8387, 147.0645, 147.5806, 147.5806, 148.3871, 149.6452, 150.3548, 151.7419, 150.8387, 151.7419, 152.6452, 153.3871, 153.1290, 151.7097, 151.7742, 152.4194, 154.1935, 152.3548, 152.6129, 152.8710, 155.7742, 158.2903, 157.9677, 153.4516, 152.7419, 152.8387, 154.3548, 155.2258, 153.6452, 154.2581, 155.8065, 153.8710, 157.7742, 156.2258, 153.1935, 152.8710, 154.0645, 155.5484, 153.2258, 152.1935, 152.5484, 152.2258, 151.5484, 152.3548, 152.5484, 154.4194, 154.9677, 153.6129, 152.5484, 152.8387, 153.2581, 154.1290, 154.9032, 152.8710, 154.2903, 152.5484, 152.6452, 153.3548, 153.1290, 151.5161, 155, 152.4839, 152.9032, 154.4516, 153.9677, 154.6452, 155.6452, 157.8710, 162.1935, 160.1935, 161.4516, 163.1613, 167.1290, 171.5806, 174.1290, 178.6774, 184.2903, 190.9032, 198.2903, 207.4839, 216.6129, 228.1290, 240.7097, 252.9032, 269.7097, 288.1290, 304.9355, 323.6774, 346.1290, 363.0645, 388.1290, 415.5484, 440.3548, 471.8710, 505.1290, 531.2581, 569.8387, 609.4839, 641.7097, 683.2903, 726.7742, 759.8065, 802.2903, 847.0645, 879.0645, 922.6774, 970.5484, 1014.7419, 1082.0968, 1163.5806, 1230.6774, 1320.1613, 1397.9677, 1430.1613, 1422.9355, 1354.3548, 1268.5484, 1153.5806, 1042, 963.5161, 872.1290, 788.8387, 731.5806, 665.2258, 607.9032, 555.9677, 496.3871, 447.2903, 411.4516, 371.2581, 335, 312.2581, 283.6129, 260.1935, 246, 224.9032, 207.9032, 198.0323, 185.4194, 175.5161, 167.9677, 160.1290, 153.8065, 148.4839, 143, 137.7097, 134.2581, 130.9032, 129.2258, 127.8710, 125.5161, 124.5806, 124.5806, 123.7742, 122.4194, 120.5806, 119.8065, 119.4516, 119.2581, 119.0645, 119.8710, 119.9355, 119.6452, 119.6774, 120.3548, 119.9032, 120.7097, 121.1613, 120.9355, 121.7419, 122.0968, 120.9677, 121.7742, 122.5161, 122.6452, 122.9355, 122.4839, 122.4194, 121.8710, 121.6774, 121.0645, 121.6774, 122.0645, 122.2903, 121.9032, 121.9032, 121.9677, 122.0323, 122.2903, 121.7742, 122.5161, 124.2581, 122.5161, 122.2581, 123.2581, 123.4516, 124.5806, 123.7742, 123.5806, 125.1290, 124.9032, 124.7097, 124.5484, 125.5484, 126.4516, 127.5161, 126.5806, 127.2903, 129.2258, 132.3226, 134.3871, 134.6452, 132.6129, 132.4516, 135.4194, 136.0645, 131, 130.9355, 131.5484, 131.9032, 132.0323, 133.0645, 133.7742, 133.6774, 134.1935, 135, 136.1935, 138.4516, 141.0323, 139.9355, 139.9355, 141.1290, 142.7419, 144.3226, 145.6452, 148.0323, 148.2258, 148.8387, 151.0968, 153.9355, 155.4516, 158.0645, 160.1613, 165.6774, 167.4839, 169.4516, 173.3548, 178.8710, 184.1290, 188.9355, 191.4839, 195.3871, 201.0645, 203.7419, 204.3226, 204.6452, 205.0323, 205.6774, 207.2258, 209.1613, 210, 212.7419, 213.6774, 213.2581, 214.7742, 216.6774, 218.2258, 220.1935, 221, 221.4194, 223.7419, 226.3871, 229.3226, 231.6129, 235.0645, 237.7097, 242.4194, 248.7097, 254.8065, 265.2258, 274.8387, 283.7097, 293.6774, 304.4839, 311.8387, 315.8710, 316.8387, 314.9677, 307.3871, 297.7097, 289.6129, 275.2903, 261.7419, 251.6452, 237.4194, 223.5806, 215.5484, 204.8387, 197.2581, 191.5161, 184.5806, 179, 175.2581, 171.5806, 168.7419, 166.4839, 164.8387, 164.0968, 164.8387, 161.4516, 160.6129, 160.3548, 159.1290, 160.3871, 160.5806, 160.5806, 159.4839, 160.3548, 160.7097, 159.3548, 159.2903, 160.0968, 160.6452, 160.7419, 160.0968, 158.8387, 158.6774, 157.6774, 157.2581, 156.8710, 157.3548, 156.3548, 156.5484, 156.4516, 155.9677, 155.4516, 155.9032, 156.1613, 156.6452, 156.1935, 156.7097, 157.7097, 159.2258, 159.9677, 161.8065, 165.9355, 159.8065, 159.6452, 159, 158.3226, 158.9032, 159, 158.0645, 157.4516, 156.7097, 156.5161, 156.2258, 157.1935, 158.7419, 158.2581, 161.3871, 161.2258, 161.0323, 159.5484, 159.2258, 159.7742, 159.0968, 158.3548, 158.4516, 158.3226, 158, 158, 157.8710, 158.2581, 157.3871, 158.9032, 158.6774, 159.9677, 160.9032, 161.1290, 161.4839, 162.6129, 161.8710, 162.3871, 162.7742, 163.3871, 164.0645, 164.5806, 164.8710, 163.1613, 164.1935, 163.9032, 163.9677, 164.5484, 164.6129, 165.2903, 165.1613, 164.7742, 164.8387, 164.7419, 166.4516, 165.0645, 165.4839, 163.4839, 165.4839, 165.3226, 164.7419, 166.8065, 167.3548, 167.5806, 167.6129, 168.4516, 169.2581, 170.0968, 170.4839, 168.2581, 168.3871, 169, 169.1290, 168.8710, 169.4839, 169.8387, 169.2581, 169.1290, 170.5484, 169.9355, 170.4839, 170.2258, 170.9677, 174.7742, 176.0323, 174.9677, 176.5484, 175.9355, 175.4194, 175.1613, 175.4516, 175.8387, 177.5484, 171.0323, 170.6129, 170.7097, 170.9032, 171.1935, 171.8065, 171.5161, 171.5806, 171.6452, 170.5484, 170.8387, 170.9032, 169.4516, 169.7742, 170.1290, 170.5161, 169.3548, 168.4516, 168.3548, 168.2581, 168.6774, 169.7419, 169.2581, 169.9032, 170.9355, 171.5806, 171.4194, 172.2581, 173.1290, 172.7097, 173.0968, 171.8710, 172.2903, 171.9032, 170.6774, 168.2581, 167.8065, 167.8387, 168.0645, 168.7742, 169.3548, 169.0323, 168.8387, 169.0968, 170.9032, 172.3226, 172.6774, 170.5484, 170.4194, 171.9355, 170.2258, 172.1290, 170.0645, 171.2903, 169.5161, 171.5806, 171.7419, 170, 169.8065, 169.7742, 170.9355, 170.6452, 176.3548, 170.5806, 170.7097, 170.2903, 170.4516, 170.8710, 171.9677, 171.2258, 170.8710, 170.1290, 170.5161, 171.0645, 174.6452, 170.6774, 169.1613, 168.7419, 166.4194, 165.9355, 167.5161, 168.9032, 165.3226, 164.2258, 163.9677, 164.6774, 166.6774, 165.8710, 165.2903, 162.8387, 162.0645, 161.9355, 162.1613, 161.9677, 161, 161.2581, 160.7419, 160.0323, 157.2581, 157.8710, 156, 154.3548, 153.9677, 153.9677, 150.8710, 149.2903, 148.6452, 146.5806, 146.8387, 146.7742, 143.9032, 141.8387, 141.8387, 142.0323, 143.8710, 146.7097, 147.4516, 147.7742, 150.7419, 155.2581, 162.9032, 170.0323, 182.3548, 200.3871, 216.8710, 244.0968, 280.2903, 314.1935, 367.2581, 434.5806, 495.0323, 579.6129, 683.1290, 765.5161, 874.7097, 987.6452, 1070.0968, 1176.1613, 1275.8065, 1337.2258, 1388.0968, 1386.7742, 1342.7742, 1243.1290, 1089.5484, 944.6774, 757.4194, 580.2258, 471.1613, 359.4839, 279.4194, 237.9677, 200.2903, 174.8710, 158.6129, 144.8710, 135.7097, 130.0323, 126.7419, 124.8065, 123.9355, 124.2581, 123.8387, 124.1935, 124.3226, 128, 127.6129, 126.0645, 126.0645, 127.1935, 127.7097, 129.3548, 129.9355, 130.1935, 130.3226, 131.4194, 132.7742, 133.6452, 134.1613, 135.7097, 138.4516, 138.3548, 137.5484, 139.0323, 139.6452, 142.4839, 143.6129, 141.5484, 143.6452, 143.5484, 144.3548, 145.7419, 144, 144.0968, 144.3871, 143.7419, 143.2903, 142.9677, 142.8710, 142.4516, 144.1290, 147.2581, 145.1290, 144.5806, 145.8387, 146.1935, 147.7097, 147.9032, 147.3548, 146.4839, 147, 145.9677, 146.0968, 146.3548, 148, 151.0323, 146.6452, 146.3226, 147.4516, 145.8387, 143.7097, 144.2258, 146.5161, 145.7097, 144.5484, 142.5806, 142.7742, 143.4516, 143.1935, 141.3548, 141.3226, 142.1935, 142.2903, 141.3226, 141.2258, 141.9355, 141.6774, 141.4516, 141.0968, 141.1290, 139.9032, 140.0323, 139.9032, 140.2903, 140.3871, 140.1613, 141.2903, 141.7097, 141.3226, 141.7419, 141.9355, 141.8710, 142.3226, 142.8065, 143.1290, 144.5806, 145.9355, 144.7742, 144.1935, 144.9032, 145.1290, 145.8710, 146.0645, 146.1613, 148.6129, 147.0968, 147.2258, 147.2258, 147.8387, 147.9032, 148.7097, 149.8065, 150.2258, 150.8065, 150.3548, 152.1290, 150.9677, 150.2903, 149.1935, 148.4194, 151.1290, 153.4516, 152.1935, 150.9677, 151.5161, 155.6452, 153.8065, 152.9677, 155.1290, 157.5806, 155.4194, 156.1935, 155.7742, 157.9355, 156.7419, 159.0968, 160.6129, 159.5484, 159.8710, 161.7419, 164.3871, 165.3871, 165.5161, 166.1935, 165.6129, 167.4839, 171.1935, 173.2581, 172.8387, 174.1290, 175.8387, 177.3548, 179.0645, 181.0968, 183.3548, 184.0323, 184.6452, 183.7097, 187.2903, 257.7097, 405.3871, 772.4194, 674, 184.2903, 170.1613, 173.5806, 176.1613, 181.0323, 184.0645, 189.8710, 195.7097, 206.7419, 211.2258, 206.2903, 204.0968, 201.3871, 198.1935, 196.7742, 195.4194, 185.0323, 177.4194, 166.9032, 157.8065, 157.9032, 157.6774, 138.4194, 122.0645);
    //yValTestOnes = newArray(380);
    //yValTestZeroes = newArray(1506-380);
    //Array.fill(yValTestOnes, 1);
    //Array.fill(yValTestZeroes, 0);
    //yVals = Array.concat(yValTestOnes, yValTestZeroes);
    //Array.print(yVals);
    yVals = getProfile();
    for (i = 0; i < yVals.length; i++) {
        yVals[i] = yVals[i] - baselineY;
    }
    binSums = newArray(binPxValues.length-1);
    for (i = 0; i < binPxValues.length-1; i ++) {
        binSums[i] = sumSingleBin(yVals, binPxValues[i], binPxValues[i+1]);
    }
    Array.print(binSums);

}

function sumSingleBin(yValArray, binLowerBound, binUpperBound) {
    slicedArray = Array.slice(yValArray, binLowerBound, binUpperBound);
    //Array.print(slicedArray);
    Array.getStatistics(slicedArray, min, max, mean, stdDev);
    return slicedArray.length * mean
}



function getBackgroundConc() {
    // laneLength = 1400;
    // selectingBaseline = false;
    // Can write this in later, should focus on the actual integration right now
    // if(selectingBaseline) {
    //     addHorizontalLineToCursor();
    // }
    // function addHorizontalLineToCursor() {
    //     while(true) {
    //         getCursorLoc(cursorX, cursorY, cursorZ, cursorModifiers);
    //         Overlay.remove;
    //         Overlay.drawLine(0, cursorY, laneLength, cursorY);
    //         Overlay.add;} 
    setTool("multi point");
    waitForUser("Select point that reflects baseline y-value of the LDL range, then press OK.");
    getSelectionCoordinates(baselineX, baselineY);
    // for (i = 0, i < baselineY.length, ) {

    // }
    getDimensions(width, height, c, z, t);
    print(height);
    print("Baseline Y:");
    Array.print(baselineY);
    print(height - baselineY[baselineY.length-1]);
    return height - baselineY[baselineY.length-1];
    }
