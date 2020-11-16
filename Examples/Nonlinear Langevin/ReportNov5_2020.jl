using Weave

using PyPlot
using Random

include("DataGen.jl") # This has many packages in it's preamble
include("../../Tools/Model_Reduction_Dev.jl")

filename = normpath(Weave.EXAMPLE_FOLDER, "ReportNov5_2020.jmd")
weave(filename; doctype = "md2pdf", out_path = :pwd)
