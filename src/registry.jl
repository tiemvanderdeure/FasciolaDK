function load_registry()
    dyr = CSV.read(joinpath(datadir, "dyr.csv"), DataFrame)
    dropmissing!(dyr, :FOEDSELSDATO) # just one row without FOEDSELSDATO
    dyrtilbes = CSV.read(joinpath(datadir, "dyrtilbes.csv"), DataFrame, types = Dict(:DATOFRA => Date, :DATOTIL => Date))
    slagt = CSV.read(joinpath(datadir, "slagtning.csv"), DataFrame)
    fund = CSV.read(joinpath(datadir, "slagtefund.csv"), DataFrame) 
    beskart = CSV.read(joinpath(datadir, "beskart.csv"), DataFrame)
    bes_chr = beskart[:, [:CHRNR, :BESNR, :CHR_ID, :BES_ID]]
    besbrugsart = CSV.read(joinpath(datadir, "besbrugsart.csv"), DataFrame)
    udstationeringer = CSV.read(joinpath(datadir, "udstationeringer.csv"), DataFrame)
#    besbrugsart_øko = besbrugsart[findall(==(16), besbrugsart.BRUGSART_ID), :]

    # recode brugsarter
    besbrugsart.øko = in.(besbrugsart.BRUGSART_ID, Ref((13,15,16,20)))
    besbrugsart.mælk = in.(besbrugsart.BRUGSART_ID, Ref((14,16,19,20)))
    # drop all burgsart > 20? These are hobbydyr, avlsdyr, etc.

    # warns but seems to be okay?
    koord = CSV.read(joinpath(datadir, "giskoord.csv"), DataFrame, decimal = ',', types = Dict(:LAT => Float64, :LON => Float64))
    # A couple of missing lon-lat values - these have XKOOR and YKOOR of 0
    dropmissing!(koord)

    return (; dyr, dyrtilbes, slagt, fund, beskart, bes_chr, besbrugsart, koord, udstationeringer)
end


