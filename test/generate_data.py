import ROOT
import numpy as np

x, y, z = np.meshgrid(np.linspace(-5, 5, 10),
                      np.linspace(-10, 10, 15),
                      np.linspace(-15, 15, 10), indexing='ij')
grid = np.transpose(np.array([x.flatten(), y.flatten(), z.flatten()]))

file = ROOT.TFile("data.root", "RECREATE")

phi = ROOT.TNtuple("phi", "phi", "x:y:z:phi")
for x, y, z in grid:
    phi.Fill(x, y, z, x + y + z)
phi.Write()

phi1 = ROOT.TNtuple("phi1", "phi1", "x:y:z:re:im")
for x, y, z in grid:
    phi1.Fill(x, y, z, x + y + z, np.sin(x + y + z))
phi1.Write()

file.Close()
