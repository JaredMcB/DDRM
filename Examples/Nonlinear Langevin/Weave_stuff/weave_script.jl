using Weave

pwd()
filename = normpath("Examples\\Nonlinear Langevin\\Weave_stuff\\Weave_test.jmd")
weave(filename;doctype = "md2pdf",
 out_path = "Examples\\Nonlinear Langevin\\Weave_stuff")
