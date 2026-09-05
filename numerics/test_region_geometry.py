"""Independent LP, planar and analytic controls for strict membership."""
import unittest
import numpy as np
from region_geometry import planar_network, lp_network, planar_patterns

class RegionChecks(unittest.TestCase):
    def compare(self,Ws):
        a,b = planar_network(Ws),lp_network(Ws)
        self.assertEqual(set(a['strict']),set(b['strict']))
        self.assertEqual(set(a['ordinary']),set(b['ordinary']))
        self.assertAlmostEqual(a['K'],b['K'],places=10)
        return a

    def test_review_counterexample(self):
        a = self.compare([np.eye(2),np.array([[1.,-1.],[-1.,2.]])])
        self.assertEqual(len(a['strict']),5)
        self.assertEqual(len(a['ordinary']),6)

    def test_dead_and_zero_rows(self):
        a = self.compare([np.eye(2),-np.ones((2,2)),np.eye(2)])
        self.assertEqual(len(a['strict']),0)
        self.assertEqual(a['K'],0)
        self.assertEqual(planar_patterns(np.array([[1.,0.],[0.,0.]])),set())

    def test_random_networks_and_width_two_law(self):
        rng = np.random.default_rng(98765)
        for L in range(1,5):
            for _ in range(12):
                Ws = [rng.normal(size=(2,2)) for _ in range(L)]
                a = self.compare(Ws)
                if L==2:
                    crossing = int(np.sum(Ws[1][:,0]*Ws[1][:,1]<0))
                    self.assertEqual(len(a['strict']),3+crossing)

    def test_narrow_cell_missed_by_uniform_grid(self):
        A = np.array([[1.,0.],[1.,1e-6]])
        self.assertEqual(len(planar_patterns(A)),4)
        t = np.linspace(0,2*np.pi,1000,endpoint=False)+.0002
        signs = A@np.array([np.cos(t),np.sin(t)])>0
        self.assertEqual(len(set(map(tuple,signs.T))),2)

if __name__=='__main__':
    unittest.main()
