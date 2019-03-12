# test import_ps function
include("..\\src\\CRISP_network.jl");
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP_Jv1.1\\data\\case6ww\\");
export_ps(ps,"data\\saved_ps\\test")
