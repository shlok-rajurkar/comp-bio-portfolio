-- Divides a single spectrum into Gaussian curves based on bins,
-- which can be edited below. By default, uses IM bins

-- Define bins based on IM parameters

local peakFolder = "/path/to/your/peaks/folder/"



local standardIMBins = {
    0, 7.65, 10.5, 14.5, 18, 19, 19.9, 
        20.49, 20.82, 21.41, 22, 22.46, 
        23.33, 25, 29.6, 33.5, 42.4, 52.29136658
}

-- Manually defined bins based on MZ correlation data

local mzCorrBins = {
    7.65, 9, 10.5, 12.5, 14.5, 15, 16, 18, 20, 23.33, 29.59

}

-- Modified IM bins for single 1200 pt sample

local modifiedIMBins = {
    3.65, 7.65, 10, 13, 18, 22, 26, 32
}

local hDLBins = {
    7.65, 9, 10.5, 13
}

local lDLBins = {

}

-- Function to clear workspace

function clearWorkspace()
    F:execute("reset")
end

-- Function to restrict range to relevant x values
function restrictRange(bins)
    F:execute("@*: A = a and not (1 < x)")
    F:execute("@*: A = a or (" .. bins[1] .. "< x and x < " .. bins[#bins] ..")")
end

-- Function to place single gaussian based on start and end values

function placeSingleCurve(start, stop, curveType)
    F:execute("@*: guess " .. curveType .. " [" .. start .. ":" .. stop .. "]")
    print(start)
    print(stop)
end

-- Fixes curve centers in place

function fixCentersInPlace(varCount)
    for i = 2, varCount, 3
    do
        F:execute("$_" .. i .."= {$_" .. i .. "}")
    end
end

-- Runs fityk command "fit" 2 times (typically enough for the fit to not change much)

function fitCurves()
    F:execute("@*: fit")
    F:execute("@*: fit")
end

-- Function to place many gaussian based on given bin values
-- Check array format

function binGaussian(bins)
    for i = 1, #bins - 1, 1
    do 
        placeSingleCurve(bins[i], bins[i+1], "Gaussian")
    end
end

function binGaussianSimple()
    for i = 1, 6, 1
    do
        F:execute("@*: guess Gaussian")
    end
end

-- Function to place many lorentzian based on given bin values

function binLorentzian(bins)
    for i = 1, #bins - 1, 1
    do 
        placeSingleCurve(bins[i], bins[i+1], "Lorentzian")
    end
end

-- Function to place many voigt curves based on given bin values

function binVoigt(bins)
    for i = 1, #bins - 1, 1
    do
        placeSingleCurve(bins[i], bins[i+1], "Voigt")
    end
end

-- Function to open CAP curve subset files

function openCAPCurveSubset(datasetCount)
    for i = 2, datasetCount + 1, 1
    do 
        F:execute("@+ < 'Z:/User Folders/SRajurkar/X3726/CAP Baseline subset n = 200/1200 pt data.csv:1:" .. i .. "::'")
    end
end



-- Optional line to export data in .peaks format
-- (tab separated values)

function exportPeaks()
    F:execute("@*: info peaks > " .. peakFolder .. "output.peaks") 
end

function fixPeaks()
    F:execute("@*: fix %.x0")
end

function analyzeCAPCurveSubsetn200()
    clearWorkspace()
    openCAPCurveSubset(200)
    restrictRange(hDLBins)
    binGaussian() -- change binning to work at center of range 
    fixCentersInPlace(200 * #hDLBins * 3)
    fitCurves()
end

-- clearWorkspace()
-- openCAPCurveSubset()
-- restrictRange(hDLBins)
-- binGaussian(hDLBins)
-- -- binGaussianSimple()
-- fitCurves()

-- -- binGaussian(mzCorrBins)
-- fitCurves()
