import numpy as np

class Run:
  """ Encapsulates output of a tinyPDR run """

  def __init__(self, filename):
    with open(filename, 'r') as f:
      header = f.readline().split()
      self.ngrid, self.ntime, self.nmols, self.nbins = map(int, header)
      self.data = np.loadtxt(f)
      self.offset_mol = 4
      self.offset_tau = self.offset_mol + self.nmols
      self.data.shape=(self.ntime, self.ngrid, self.offset_mol+self.nmols+self.nbins)

  def getSnapshot(self, itime):
    return self.data[itime,:,:]

  def getMolecules(self, itime=None):
    mols = self.data[:,:,self.offset_mol:self.offset_mol+self.nmols]
    if itime:
      return mols[itime,:,:]
    return mols

  def getTau(self, itime=None):
    tau = self.data[:,:,self.offset_tau:]
    if itime:
      return tau[itime,:,:]
    return tau

  def getTime(self, itime=None):
    time = self.data[:,:,0]
    if itime:
      return time[itime,:]
    return time

  def getR(self, itime=None):
    r = self.data[:,:,1]
    if itime:
      return r[itime,:]
    return r

  def getAv(self, itime=None):
    Av = self.data[:,:,2]
    if itime:
      return Av[itime,:]
    return Av

  def getTgas(self, itime=None):
    Tgas = self.data[:,:,3]
    if itime:
      return Tgas[itime,:]
    return Tgas