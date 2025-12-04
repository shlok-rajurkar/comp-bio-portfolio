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
    7.65, 10, 13, 18, 22, 32
}

-- Function to restrict range to relevant x values
function restrictRange()
    F:execute("@*: A = a and not (1 < x)")
    F:execute("@*: A = a or (7.65 < x and x < 32)")
end

-- Function to place single gaussian based on start and end values

function placeSingleCurve(start, stop, curveType)
    F:execute("@*: guess " .. curveType .. " [" .. start .. ":" .. stop .. "]")
    print(start)
    print(stop)
end

-- Runs fityk command "fit" 5 times (typically enough for the fit to not change much)

function fitCurves()
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

-- Function to place many lorentzian based on given bin values

function binLorentzian(bins)
    for i = 1, #bins - 1, 1
    do 
        placeSingleCurve(bins[i], bins[i+1], "Lorentzian")
    end
end

-- Function to open CAP curve subset files

function openCAPCurveSubset()
    for i = 2, 201, 1
    do 
        F:execute("@+ < 'Z:/User Folders/SRajurkar/X3726/CAP Baseline subset n = 200/1200 pt data.csv:1:" .. i .. "::'")
    end
end



-- Optional line to export data in .peaks format
-- (tab separated values)

function exportPeaks()
    F:execute("@*: info peaks > " .. peakFolder .. "output.peaks") 
end



openCAPCurveSubset()
restrictRange()
binLorentzian(modifiedIMBins)
fitCurves()
-- binGaussian(mzCorrBins)
-- fitCurves()
