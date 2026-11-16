using JutulDarcy.GeoEnergyIO: BinaryKeyword, write_binary
using JutulDarcy.GeoEnergyIO.StandardIO: byteorder, keyword_type, read_record, write_record
using JutulDarcy.StandardOutput: grid_layout

function restart_well_records(case)
    case.model isa MultiModel || return nothing
    wells = [(String(name), unique(well.perforations.reservoir))
        for (name, well) in JutulDarcy.get_model_wells(case)
        if hasproperty(well, :perforations) && !isempty(well.perforations.reservoir)]
    isempty(wells) && return nothing

    (; dims, cell_map) = grid_layout(case)
    logical_cells = CartesianIndices(dims)
    nwells = length(wells)
    max_connections = maximum(length(cells) for (_, cells) in wells)
    # Eclipse restart widths for IWEL, ZWEL and ICON used by OPM Flow.
    iwel = zeros(Int32, 155, nwells)
    zwel = fill("", 3, nwells)
    icon = zeros(Int32, 26, max_connections, nwells)
    for (w, (name, cells)) in enumerate(wells)
        ijk = logical_cells[cell_map[cells]]
        i, j, _ = Tuple(first(ijk))
        first_k, last_k = extrema(cell -> cell[3], ijk)
        iwel[1:5, w] .= (i, j, first_k, last_k, length(ijk))
        zwel[1, w] = name
        for (connection, cell) in enumerate(ijk)
            icon[1:4, connection, w] .= (connection, Tuple(cell)...)
        end
    end

    # INTEHEAD positions for NWELLS, NCWMAX, NIWELZ, NZWELZ and NICONZ.
    dimensions = (nwells, max_connections, size(iwel, 1), size(zwel, 1), size(icon, 1))
    records = (BinaryKeyword("IWEL", iwel), BinaryKeyword("ZWEL", zwel), BinaryKeyword("ICON", icon))
    return (; dimensions, records)
end

function add_restart_wells(path, case)
    wells = restart_well_records(case)
    isnothing(wells) && return nothing
    # Work on the staged file. Stream the reservoir arrays through unchanged.
    mktemp(dirname(path)) do temporary, dest
        open(path, "r") do source
            copy_restart_with_wells(dest, source, wells)
        end
        close(dest)
        Base.Filesystem.rename(temporary, path)
    end
    return nothing
end

function copy_restart_with_wells(dest, source, wells)
    endian = :big # GeoEnergyIO's default for write_standard_output.
    add_wells = false
    while !eof(source)
        header = read_record(source, endian; expected = 16)
        name = rstrip(String(header[1:8]))
        count = Int(byteorder(only(reinterpret(Int32, header[9:12])), endian))
        if name == "INTEHEAD"
            data = read_record(source, endian; expected = count * sizeof(Int32))
            fields = reinterpret(Int32, data)
            map!(x -> byteorder(x, endian), fields, fields)
            add_wells = iszero(fields[17])
            if add_wells
                for (index, value) in zip((17, 18, 25, 28, 33), wells.dimensions)
                    fields[index] = value
                end
            end
            write_binary(dest, (BinaryKeyword(name, fields),); endian)
        else
            if add_wells && name == "STARTSOL"
                write_binary(dest, wells.records; endian)
            end
            _, width = keyword_type(String(header[13:16]))
            write_record(dest, header, endian)
            remaining = count * width
            while remaining > 0
                data = read_record(source, endian; limit = remaining)
                isempty(data) && throw(ArgumentError("Empty payload in restart keyword $name"))
                write_record(dest, data, endian)
                remaining -= length(data)
            end
        end
    end
    return nothing
end
