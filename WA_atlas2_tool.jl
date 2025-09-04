using Rasters
using ArchGDAL
using Interpolations
using DataFrames
using PrettyTables
using DimensionalData: X, Y, At, Near

using Gtk4

import StatsBase
import Interpolations: Linear, Nearest
import ArchGDAL


using Gtk4
using DataFrames
using PrettyTables
import PrettyTables: tf_ascii, ft_round


"""
"""
function generate_ascii_name(state::String, recurrence::Int, duration::Int)
    state = lowercase(state)
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
            depth=val_bilinear,
            intensity=val_bilinear/duration)

end

"""
"""
function generate_one_hour_data(six_hour_data, twentyfour_hour_data, region::Int, lat::Float64, long::Float64, elev::Float64)
        
        x_one = six_hour_data.two_year
        x_two = twentyfour_hour_data.two_year

        x_three = six_hour_data.hundred_year
        x_four = twentyfour_hour_data.hundred_year

        x_five = lat - 40
        x_six = abs(long) - 100
        zee = elev / 100
        
    if region === 1
        two_one = 0.019 + 0.711*(x_one*(x_one/x_two)) + 0.001*zee
        hundred_one = 0.338 + 0.670*(x_three*(x_three/x_four)) + 0.001*zee
    elseif region === 2
        two_one = 0.077 + 0.715*(x_one*(x_two/x_three)) - 0.0004*(x_five/x_six)
        hundred_one = 0.187 + 0.833*(x_three*(x_three/x_four))
    elseif region === 3
        two_one = 0.157 + 0.513*(x_one*(x_one/x_two))
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
    six_hours = general_recurrence_nomogram(twosix.depth, hundredsix.depth)
    println(six_hours)

    # get all 24-hour depths
    twentyfour_hours = general_recurrence_nomogram(twotwentyfour.depth, hundredtwentyfour.depth)

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
                        ("24-hr" => twentyfour_hours)]
    out_df = DataFrame(recurrence=String[], duration=String[], depth=Float64[])
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
    #pretty_results = pretty_table(wide_df)
    return wide_df
end




# testing
#test = generate_full_data("wa", 3, 47.115047, -123.754755, 341.0)


win = GtkWindow("Washington Atlas2 Vol9 Tool", 600, 600)

# create box container
vbox = GtkBox(:v)
set_gtk_property!(vbox, :spacing, 10)
set_gtk_property!(vbox, :margin_left, 20)
set_gtk_property!(vbox, :margin_right, 20)
set_gtk_property!(vbox, :margin_top, 20)
set_gtk_property!(vbox, :margin_bottom, 20)


# Create input fields
# state
state_label = GtkLabel("State:")
set_gtk_property!(state_label, :halign, Gtk4.GtkAlign.START)
state_ent = GtkEntry()
set_gtk_property!(state_ent, :text, "WA")

# region
region_label = GtkLabel("Region:")
set_gtk_property!(region_label, :halign, Gtk4.GtkAlign.START)
region_ent = GtkEntry()
set_gtk_property!(region_ent, :text, "3")

# latitude
lat_label = GtkLabel("Latitude:")
set_gtk_property!(lat_label, :halign, Gtk4.GtkAlign.START)
lat_ent = GtkEntry()
set_gtk_property!(lat_ent, :text, "47.115078")

# longitude
long_label = GtkLabel("Longitude:")
set_gtk_property!(long_label, :halign, Gtk4.GtkAlign.START)
long_ent = GtkEntry()
set_gtk_property!(long_ent, :text, "-123.754755")

# Elevation
elev_label = GtkLabel("Elevation:")
set_gtk_property!(elev_label, :halign, Gtk4.GtkAlign.START)
elev_ent = GtkEntry()
set_gtk_property!(elev_ent, :text, "341")

# Results area
results_label = GtkLabel("Results:")
set_gtk_property!(results_label, :halign, Gtk4.GtkAlign.START)
results_textview = GtkTextView()
results_buffer = GtkTextBuffer()
set_gtk_property!(results_textview, :monospace, true)
set_gtk_property!(results_textview, :buffer, results_buffer)
set_gtk_property!(results_textview, :editable, false)
set_gtk_property!(results_textview, :wrap_mode, Gtk4.GtkWrapMode.NONE)

# scrolled window for results
results_scrolled = GtkScrolledWindow()
push!(results_scrolled, results_textview)
set_gtk_property!(results_scrolled, :hscrollbar_policy, Gtk4.GtkPolicyType.ALWAYS)
set_gtk_property!(results_scrolled, :vscrollbar_policy, Gtk4.GtkPolicyType.AUTOMATIC)
set_gtk_property!(results_scrolled, :height_request, 200)

# run button
run_button = GtkButton("Compute Rainfall")

function on_run(w)
    try
        state = get_gtk_property(state_ent, :text, String)
        region = parse(Int, get_gtk_property(region_ent, :text, String))
        lat = parse(Float64, get_gtk_property(lat_ent, :text, String))
        long = parse(Float64, get_gtk_property(long_ent, :text, String))
        elev = parse(Float64, get_gtk_property(elev_ent, :text, String))

        # clear prev results
        set_gtk_property!(results_buffer, :text, "Computing...")

        # call function
        results = generate_full_data(state, region, lat, long, elev)

        # format and display
        if results !== nothing
            try
                result_text = sprint() do io
                    pretty_table(io, results;
                                backend=Val(:text),
                                alignment=:c,
                                crop=:none,
                                formatters=PrettyTables.ft_round(2),
                                tf=tf_ascii,
                                header=(names(results), nothing))
                end
                set_gtk_property!(results_buffer, :text, result_text)
            catch e
                @warn "PrettyTables formatting failed: $e"
                result_text = string(results)
                set_gtk_property!(results_buffer, :text, result_text)
            end
            
        else
            set_gtk_property!(results_buffer, :text, "Error: Could not compute")            
        end

    catch e
        error_text = "Error: $e\n\nPlease check your input values"
        set_gtk_property!(results_buffer, :text, error_text)
        @warn "Error in computation: $e"
    end
end
signal_connect(on_run, run_button, "clicked")

# Add all widgets to the vertical box
push!(vbox, state_label)
push!(vbox, state_ent)
push!(vbox, region_label) 
push!(vbox, region_ent)
push!(vbox, lat_label)
push!(vbox, lat_ent)
push!(vbox, long_label)
push!(vbox, long_ent)
push!(vbox, elev_label)
push!(vbox, elev_ent)
push!(vbox, run_button)
push!(vbox, results_label)
push!(vbox, results_scrolled)

# Add the box to the window
push!(win, vbox)

# Replace the existing signal_connect for destroy with this:
signal_connect(win, "destroy") do widget
    Gtk4.gtk_quit()
    exit()
end

showall(win)

# Keep the main thread alive
Gtk4.gtk_main()
