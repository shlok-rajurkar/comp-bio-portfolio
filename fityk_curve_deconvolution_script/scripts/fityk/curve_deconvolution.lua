-- Divides a single spectrum into Gaussian curves based on bins,
-- which can be edited below. By default, uses IM bins

-- Define bins based on IM parameters

local standardIMBins = {
    0, 7.65, 10.5, 14.5, 18, 19, 19.9, 
        20.49, 20.82, 21.41, 22, 22.46, 
        23.33, 25, 29.6, 33.5, 42.4, 52.29136658
}

-- Manually defined bins based on MZ correlation data

local mzCorrBins = {
    7.65, 9, 10.5, 12.5, 14.5, 15, 16, 18, 20, 23.33, 29.59

}

-- Need to restrict range in a single line because later calls
-- overwrite previous ones...

function restrictRange()
    F:execute("@*: A = a and not (1 < x)")
    F:execute("@*: A = a or (7.65 < x and x < 29.59)")
end

-- Function to place single gaussian based on start and end values

function placeSingleGaussian(start, stop)
    F:execute("@*: guess Gaussian [" .. start .. ":" .. stop .. "]")
    print(start)
    print(stop)
end

-- Runs fityk command "fit" 5 times (typically enough for the fit to not change much)

function fitCurves()
    for _ = 1, 5, 1
    do
        F:execute("@*: fit")
    end
end

-- Function to place many gaussian based on given bin values
-- Check array format

function binGaussian(bins)
    for i = 1, #bins - 1, 1
    do 
        placeSingleGaussian(bins[i], bins[i+1])
    end
end

restrictRange()
binGaussian(mzCorrBins)
fitCurves()

-- Optional line to export data in .peaks format
-- (tab separated values)

function exportPeaks()
    F:execute("@*: info peaks > 'Z:/User Folders/SRajurkar/X3726/Peak Data/'")
end