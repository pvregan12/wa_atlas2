using Rasters
using ArchGDAL
using Interpolations
using DataFrames
using DimensionalData: X, Y, At, Near

import StatsBase
import Interpolations: Linear, Nearest
import ArchGDAL

"""
"""
function generate_ascii_name(state::String, recurrence::Int, duration::Int)
    return "na2_"*state*"_"*string(recurrence)*"yr"*string(duration)*"hr.asc"
end

"""
"""
function read_ascii_grid(state::String, recurrence::Int, duration::Int, lat::Float64, long::Float64)
    # get ascii file name
    ascii_name = generate_ascii_name(state, recurrence, duration)
    script_dir = @__DIR__

    # initialize rast outside try/catch block
    rast = nothing
    # pull ascii file
    try
        println(joinpath(script_dir, ascii_name))
        rast = Raster(joinpath(script_dir, ascii_name))
    catch e
        @error "Issue finding ascii grid named $ascii_name : $e"
        return nothing
    end  
    # bilinear interpolation between lat/long
    #ri = interpolate(rast, BSpline(Linear()))
    #pt = (X(At(long)), Y(At(lat)))
    #val_bilinear = ri[X(At(long)), Y(At(lat))]
    #val_bilinear = Rasters.extract(rast, (X(At(long)), Y(At(lat))); method=Linear())
    
    val_bilinear = rast[X(Near(long)), Y(Near(lat))]
    
    
    if val_bilinear < 0 || isnan(val_bilinear)
        val_bilinear = NaN
    else
        # divide by 100,000 to get to inches of depth
        val_bilinear = val_bilinear / 100000
    end
    println("$recurrence year $duration hour depth: $val_bilinear")
    return (state=state,
            recurrence=recurrence,
            duration=duration,
            latitude=lat,
            longitude=long,
            depth=val_bilinear)

end

"""
"""
function generate_one_hour_data(six_hour_data, twentyfour_hour_data, region::Int, lat::Float64, long::Float64, elev::Float64)
        
        x_one = six_hour_data.two_year
        x_two = twentyfour_hour_data.two_year

        x_three = six_hour_data.hundred_year
        x_four = twentyfour_hour_data.hundred_year

        x_five = lat - 40
        x_six = long - 100
        zee = elev / 100
        
    if region === 1
        two_one = 0.019 + 0.711*(x_one*(x_one/x_two)) + 0.001*zee
        hundred_one = 0.338 + 0.670*(x_three*(x_three/x_four)) + 0.001*zee
    elseif region === 2
        two_one = 0.077 + 0.715*(x_one*(x_two/x_three)) - 0.0004*(x_five/x_six)
        hundred_one = 0.187 + 0.833*(x_three*(x_three/x_four))
    elseif region === 3
        two_one = 0.157 + 0.513*(x_one*(x_one.x_two))
        hundred_one = 0.324 + 0.752*(two_one*(x_three/x_one))
    elseif region === 4
        two_one = 0.160 + 0.520*(x_one*(x_one/x_two))
        hundred_one = 0.177 + 0.965*(two_one*(x_three/x_one))
    else
        @error "Invalid region number: $region"
    end
    one_hour = general_recurrence_nomogram(two_one, hundred_one)
    return one_hour
end

"""
"""
function compute_twelve_hour_depths(six_hours, twentyfour_hours)
    output_tuple = (
                    two_year = (six_hours.two_year + twentyfour_hours.two_year)/2,
                    five_year = (six_hours.five_year + twentyfour_hours.five_year)/2,
                    ten_year = (six_hours.ten_year + twentyfour_hours.ten_year)/2,
                    twentyfive_year = (six_hours.twentyfive_year + twentyfour_hours.twentyfive_year)/2,
                    fifty_year = (six_hours.fifty_year + twentyfour_hours.fifty_year)/2,
                    hundred_year = (six_hours.hundred_year + twentyfour_hours.hundred_year)/2
    )
    return output_tuple
end

"""
"""
function compute_two_hour_depths(one_hours, six_hours, region::Int)
    if region === 1
        return (two_year = 0.250*six_hours.two_year + 0.750*one_hours.two_year,
            five_year = 0.250*six_hours.five_year + 0.750*one_hours.five_year,
            ten_year = 0.250*six_hours.ten_year + 0.750*one_hours.ten_year,
            twentyfive_year = 0.250*six_hours.twentyfive_year + 0.750*one_hours.twentyfive_year,
            fifty_year = 0.250*six_hours.fifty_year + 0.750*one_hours.fifty_year,
            hundred_year = 0.250*six_hours.hundred_year + 0.750*one_hours.hundred_year
       )
    elseif region === 2
       return (two_year = 0.278*six_hours.two_year + 0.722*one_hours.two_year,
            five_year = 0.278*six_hours.five_year + 0.722*one_hours.five_year,
            ten_year = 0.278*six_hours.ten_year + 0.722*one_hours.ten_year,
            twentyfive_year = 0.278*six_hours.twentyfive_year + 0.722*one_hours.twentyfive_year,
            fifty_year = 0.278*six_hours.fifty_year + 0.722*one_hours.fifty_year,
            hundred_year = 0.278*six_hours.hundred_year + 0.722*one_hours.hundred_year
       )
    elseif region === 3 || region === 4
       return (two_year = 0.240*six_hours.two_year + 0.760*one_hours.two_year,
            five_year = 0.240*six_hours.five_year + 0.760*one_hours.five_year,
            ten_year = 0.240*six_hours.ten_year + 0.760*one_hours.ten_year,
            twentyfive_year = 0.240*six_hours.twentyfive_year + 0.760*one_hours.twentyfive_year,
            fifty_year = 0.240*six_hours.fifty_year + 0.760*one_hours.fifty_year,
            hundred_year = 0.240*six_hours.hundred_year + 0.760*one_hours.hundred_year
       )
    else
        @error "Invalid region number: $region"
    end

end

"""
"""
function compute_three_hour_depths(one_hours, six_hours, region::Int)
    if region === 1
        return (two_year = 0.467*six_hours.two_year + 0.533*one_hours.two_year,
            five_year = 0.467*six_hours.five_year + 0.533*one_hours.five_year,
            ten_year = 0.467*six_hours.ten_year + 0.533*one_hours.ten_year,
            twentyfive_year = 0.467*six_hours.twentyfive_year + 0.533*one_hours.twentyfive_year,
            fifty_year = 0.467*six_hours.fifty_year + 0.533*one_hours.fifty_year,
            hundred_year = 0.467*six_hours.hundred_year + 0.533*one_hours.hundred_year
       )
    elseif region === 2
       return (two_year = 0.503*six_hours.two_year + 0.497*one_hours.two_year,
            five_year = 0.503*six_hours.five_year + 0.497*one_hours.five_year,
            ten_year = 0.503*six_hours.ten_year + 0.497*one_hours.ten_year,
            twentyfive_year = 0.503*six_hours.twentyfive_year + 0.497*one_hours.twentyfive_year,
            fifty_year = 0.503*six_hours.fifty_year + 0.497*one_hours.fifty_year,
            hundred_year = 0.503*six_hours.hundred_year + 0.497*one_hours.hundred_year
       )
    elseif region === 3 || region === 4
       return (two_year = 0.468*six_hours.two_year + 0.532*one_hours.two_year,
            five_year = 0.468*six_hours.five_year + 0.532*one_hours.five_year,
            ten_year = 0.468*six_hours.ten_year + 0.532*one_hours.ten_year,
            twentyfive_year = 0.468*six_hours.twentyfive_year + 0.532*one_hours.twentyfive_year,
            fifty_year = 0.468*six_hours.fifty_year + 0.532*one_hours.fifty_year,
            hundred_year = 0.468*six_hours.hundred_year + 0.532*one_hours.hundred_year
       )
    else
        @error "Invalid region number: $region"
    end

end

"""
"""
function general_recurrence_nomogram(two_year, hundred_year)
    mapping_array = [1 0.75073 0.58722 0.37898 0.18518 0;
               0 0.24927 0.41278 0.62102 0.81482 1]

    output_tuple = (
                    two_year = two_year,
                    five_year = (two_year*mapping_array[1,2]) + (hundred_year*mapping_array[2,2]),
                    ten_year = (two_year*mapping_array[1,3]) + (hundred_year*mapping_array[2,3]),
                    twentyfive_year = (two_year*mapping_array[1,4]) + (hundred_year*mapping_array[2,4]),
                    fifty_year = (two_year*mapping_array[1,5]) + (hundred_year*mapping_array[2,5]),
                    hundred_year = hundred_year
    )

    return output_tuple
end

"""
"""
function generate_full_data(state::String, region::Int, lat::Float64, long::Float64, elev::Float64)

    # first, get 2-year and 100-year 6-hour and 24-hour depths
    twosix = read_ascii_grid(state, 2, 6, lat, long)
    twotwentyfour = read_ascii_grid(state, 2, 24, lat, long)
    hundredsix = read_ascii_grid(state, 100, 6, lat, long)
    hundredtwentyfour = read_ascii_grid(state, 100, 24, lat, long)

    # get all 6-hour depths, including storing 2 and 24 hour
    six_hours = general_recurrence_nomogram(twosix, hundredsix)

    # get all 24-hour depths
    twentyfour_hours = general_recurrence_nomogram(twotwentyfour, hundredtwentyfour)

    # get 1-hour depths
    one_hours = generate_one_hour_data(six_hours, twentyfour_hours, region, lat, long, elev)

    # compute the 2 and 3 hour depths
    two_hours = compute_two_hour_depths(one_hours, six_hours, region)
    three_hours = compute_three_hour_depths(one_hours, six_hours, region)

    # get all 12-hour depths
    twelve_hours = compute_twelve_hour_depths(six_hours, twentyfour_hours)

    duration_labels = [("1-hr" => one_hours),
                        ("2-hr" => two_hours),
                        ("3-hr" => three_hours),
                        ("6-hr" => six_hours),
                        ("12-hr" => twelve_hours),
                        ("24-hr" => twentyfour_hours_hours)]
    out_df = DataFrame(recurrence=Int[], duration=String[], depth=Float64[])
    # format into a dataframe
    for (dur_label, nt) in duration_labels
        for (k, v) in pairs(nt)
            rec = String(k)
            depth = v isa AbstractVector ? only(v) : Float64(v)
            push!(out_df, (rec, dur_label, depth))
        end
    end
    # long => wide format
    wide_df = unstack(out_df, :recurrence, :duration, :depth)
    return wide_df
end




# testing
test = generate_full_data("wa", 3, 47.115047, -123.754755, 341.0)


