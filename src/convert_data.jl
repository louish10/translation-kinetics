using CSV, DataFrames

df = DataFrame(CSV.File("data/ExperimentalData.csv"))


select!(
    df,
    "Protein IDs" => "id",
    "Protein copy number average [molecules/cell]" => "pNumber",
    "mRNA copy number average [molecules/cell]" => "mNumber",
    "transcription rate (vsr) average [molecules/(cell*h)]" => "alpha",
    "translation rate constant (ksp) average [molecules/(mRNA*h)]" => "beta",
)

replace!(df.pNumber, "NA" => "")
replace!(df.mNumber, "NA" => "")
replace!(df.alpha, "NA" => "")
replace!(df.beta, "NA" => "")
#df["alpha"] = parse.([Float64], df[!, "alpha"])
#df["beta"] = parse.([Float64], df[!, "beta"])

CSV.write("data/experimental_data_sanitised.csv", df)
