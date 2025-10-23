-- Divides a single spectrum into Gaussian curves based on bins,
-- which can be edited below. By default, uses 

-- Define bins based on IM parameters

local bins = {0, 7.65, 10.5, 14.5, 18, 19, 19.9,
        20.49, 20.82, 21.41, 22, 22.46, 
        23.33, 25, 29.6, 33.5, 42.4, 54.7, 100}

-- Function to place single gaussian based on start and end values

function placeSingleGaussian(start, stop)
    F:execute("@*: guess Gaussian [" .. start .. ":" .. stop .. "]")
    print(start)
    print(stop)
end

