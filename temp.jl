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
                                alignment=:c,
                                crop=:none,
                                formatters=ft_round(4))
                end
            catch e
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