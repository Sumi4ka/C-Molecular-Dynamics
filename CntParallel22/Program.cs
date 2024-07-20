using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Threading;

namespace CntCnt
{
    abstract class AtomRungeKutta4Coof
    {
        public double x0, y0, z0;
        public double xb1, yb1, zb1;
        public double xb2, yb2, zb2;
        public double xb3, yb3, zb3;
        public double xb4, yb4, zb4;
    }
    abstract class Atom : AtomRungeKutta4Coof
    {
        readonly Molecule M;
        public double x, y, z, u, v, w, m, eps = 0, sigma = 0, xA, yA, zA, ro, eps0, sigma0, c, U, FLJ;
        public Atom(double x, double y, double z, Molecule M) { this.x = x; this.y = y; this.z = z; this.M = M; }
        public double FLJx, FLJy, FLJz;
        public void VelocityStep()
        {
            u = M.OmegaY * z - M.OmegaZ * y;
            v = M.OmegaZ * x - M.OmegaX * z;
            w = M.OmegaX * y - M.OmegaY * x;
        }
        public bool flag = false;
        public void ParallelMethod()
        {
            FLJx = 0.0; FLJy = 0.0; FLJz = 0.0;
            foreach (Molecule m in M.ExternalObject)
            {
                foreach (Atom a1 in m.Atoms)
                {
                    ro = Math.Sqrt(Math.Pow(xA - a1.xA, 2.0) + Math.Pow(yA - a1.yA, 2.0) + Math.Pow(zA - a1.zA, 2.0));
                    eps0 = Math.Sqrt(this.eps * a1.eps) * Constant.KB;
                    sigma0 = (this.sigma + a1.sigma) / 2.0;
                    c = Math.Pow(sigma0 / ro, 6.0);
                    FLJ = 24.0 * eps0 / ro / ro * c * (2.0 * c - 1.0);
                    FLJx += FLJ * (xA - a1.xA);
                    FLJy += FLJ * (yA - a1.yA);
                    FLJz += FLJ * (zA - a1.zA);
                }
            }
            flag = true;
        }
        public void ParallelMethod2()
        {
            FLJx = 0.0; FLJy = 0.0; FLJz = 0.0;
            foreach (Molecule m in M.ExternalObject)
            {
                foreach (Atom a1 in m.Atoms)
                {
                    ro = Math.Sqrt(Math.Pow(xA - a1.xA, 2.0) + Math.Pow(yA - a1.yA, 2.0) + Math.Pow(zA - a1.zA, 2.0));
                    eps0 = Math.Sqrt(this.eps * a1.eps) * Constant.KB;
                    sigma0 = (this.sigma + a1.sigma) / 2.0;
                    c = Math.Pow(sigma0 / ro, 6.0);
                    U = 4.0 * eps0 * c * (c - 1.0);
                }
            }
            flag = true;
        }
    }
    sealed class Carbon : Atom
    {
        public Carbon(double x, double y, double z, Molecule M) : base(x, y, z, M)
        {
            m = 12.011 * Constant.atomMass;
            eps = 12.5;
            sigma = 0.34;
        }
    }
    abstract class MoleculeRungeKutta4Coof
    {
        public double u0, v0, w0, x0, y0, z0, Kx0, Ky0, Kz0;
        public double ub1, vb1, wb1, xb1, yb1, zb1, Lx1, Ly1, Lz1;
        public double ub2, vb2, wb2, xb2, yb2, zb2, Lx2, Ly2, Lz2;
        public double ub3, vb3, wb3, xb3, yb3, zb3, Lx3, Ly3, Lz3;
        public double ub4, vb4, wb4, xb4, yb4, zb4, Lx4, Ly4, Lz4;
    }
    sealed class Molecule : MoleculeRungeKutta4Coof
    {
        public double x, y, z, u, v, w, OmegaX, OmegaY, OmegaZ, Mass = 0.0;
        public List<Atom> Atoms = new List<Atom>();
        public string name;
        public Molecule(string name) { this.name = name; }
        public void AddAtom(Atom a) { Atoms.Add(a); this.Mass += a.m; }
        public void InitiateSpaceCoor(double x, double y, double z) { this.xM0 = this.x = x; this.yM0 = this.y = y; this.zM0 = this.z = z; }
        public void InitiateVelocityCoor(double u, double v, double w) { this.uM0 = this.u = u; this.vM0 = this.v = v; this.wM0 = this.w = w; }
        public void InitiateAngleVelocityCoor(double OmegaX, double OmegaY, double OmegaZ) { this.OmegaX0 = this.OmegaX = OmegaX; this.OmegaY0 = this.OmegaY = OmegaY; this.OmegaZ0 = this.OmegaZ = OmegaZ; }
        //External
        public readonly List<Molecule> ExternalObject = new List<Molecule>();
        public void AddEObject(Molecule m) { ExternalObject.Add(m); }
        //Supporting coofs and methods
        public double A, B, C, D, E, F, Kx, Ky, Kz, ub, vb, wb, Lx, Ly, Lz;
        public double xM0, yM0, zM0, uM0, vM0, wM0, OmegaX0, OmegaY0, OmegaZ0;
        public void CMToZero() //функция переноса центра масс в нуль сис.коор.
        {
            double xM = 0.0, yM = 0.0, zM = 0.0;
            foreach (Atom a in Atoms)
            {
                xM += a.m * a.x;
                yM += a.m * a.y;
                zM += a.m * a.z;
            }
            xM /= Mass;
            yM /= Mass;
            zM /= Mass;
            foreach (Atom a in Atoms) { a.x -= xM; a.y -= yM; a.z -= zM; }
        }
        public void CoofTenzor()
        {
            A = 0.0;
            B = 0.0;
            C = 0.0;
            D = 0.0;
            E = 0.0;
            F = 0.0;
            for (int j = 0; j < Atoms.Count; j++)
            {
                A += Atoms[j].m * (Math.Pow(Atoms[j].y, 2) + Math.Pow(Atoms[j].z, 2));
                B += Atoms[j].m * (Math.Pow(Atoms[j].x, 2) + Math.Pow(Atoms[j].z, 2));
                C += Atoms[j].m * (Math.Pow(Atoms[j].x, 2) + Math.Pow(Atoms[j].y, 2));

                D -= Atoms[j].m * Atoms[j].y * Atoms[j].z;
                E -= Atoms[j].m * Atoms[j].x * Atoms[j].z;
                F -= Atoms[j].m * Atoms[j].x * Atoms[j].y;
            }
        }
        public void AngleVelocity()
        {
            double delta0 = A * B * C + E * F * D + E * F * D - B * E * E - C * F * F - A * D * D;
            if (delta0 != 0.0)
            {
                double delta1 = Kx * (B * C - D * D) + Ky * (E * D - F * C) + Kz * (F * D - B * E);
                double delta2 = Kx * (E * D - F * C) + Ky * (A * C - E * E) + Kz * (F * E - A * D);
                double delta3 = Kx * (F * D - E * B) + Ky * (E * F - A * D) + Kz * (A * B - F * F);
                OmegaX = delta1 / delta0; OmegaY = delta2 / delta0; OmegaZ = delta3 / delta0;
            }
            else { OmegaX = 0.0; OmegaY = 0.0; OmegaZ = 0.0; }
        }
        public void KineticMoment()
        {
            Kx = A * OmegaX + F * OmegaY + E * OmegaZ;
            Ky = F * OmegaX + B * OmegaY + D * OmegaZ;
            Kz = E * OmegaX + D * OmegaY + C * OmegaZ;
        }
        public void MomentForces()
        {
            Lx = Ly = Lz = 0.0;
            foreach (Atom a in Atoms)
            {
                Lx += a.y * a.FLJz - a.z * a.FLJy;
                Ly += a.z * a.FLJx - a.x * a.FLJz;
                Lz += a.x * a.FLJy - a.y * a.FLJx;
            }
        }
    }
    sealed class Solver
    {
        public double dt; public int T;
        public Solver(int T, double dt)
        {
            this.dt = dt; this.T = T;
            Data.SetTdt(T, dt);
        }
        public readonly List<Molecule> Molecules = new List<Molecule>();
        public void AddMolecule(Molecule M)
        {
            foreach (Molecule m in Molecules)
            {
                m.AddEObject(M);
                M.AddEObject(m);
            }
            Molecules.Add(M);
        }
        public void UVW()
        {
            foreach (Molecule m in Molecules)
            {
                m.ub = 0.0; m.vb = 0.0; m.wb = 0.0;
            }
            foreach (Molecule m in Molecules) foreach (Atom a in m.Atoms)
                {
                    ThreadPool.QueueUserWorkItem((i) =>
                    {
                        a.ParallelMethod();
                    });
                }
            foreach (Molecule m in Molecules)
            {
                bool flag = true;
                while (flag)
                {
                    foreach (Atom a in m.Atoms)
                    {
                        flag = flag && a.flag;
                    }
                    flag = !flag;
                }
                foreach (Atom a in m.Atoms) a.flag = false;
            }
            foreach (Molecule m in Molecules)
            {
                foreach (Atom a in m.Atoms)
                {
                    m.ub += a.FLJx;
                    m.vb += a.FLJy;
                    m.wb += a.FLJz;
                }
                m.ub /= m.Mass; m.vb /= m.Mass; m.wb /= m.Mass;
            }
        }
        public double Balance()
        {
            double H = 0;
            foreach (Molecule m in Molecules)
            {
                double Omega = Math.Sqrt(Math.Pow(m.OmegaX, 2.0) + Math.Pow(m.OmegaY, 2.0) + Math.Pow(m.OmegaZ, 2.0));
                if (Omega == 0)
                {
                    double alpha = m.OmegaX/Omega, beta = m.OmegaY / Omega, gamma = m.OmegaZ / Omega;
                    double J = (m.A * Math.Pow(alpha, 2) + m.B * Math.Pow(beta, 2) + m.C * Math.Pow(gamma, 2) + 2 * m.D * beta * gamma + 2 * m.E * gamma * alpha + 2 * m.F * alpha * beta);
                    H += J * Math.Pow(Omega, 2) / 2.0;
                }
                double E = m.Mass * (Math.Pow(m.u, 2.0) + Math.Pow(m.v, 2.0) + Math.Pow(m.w, 2.0)) / 2;
                H += E;
                double U = 0.0;
                foreach (Atom a in m.Atoms)
                {
                    ThreadPool.QueueUserWorkItem((i) =>
                    {
                         a.ParallelMethod2();
                    });
                }
                bool flag = true;
                while (flag)
                {
                    foreach (Atom a in m.Atoms) flag = flag && a.flag;
                    flag = !flag;
                }
                foreach (Atom a in m.Atoms) a.flag = false;
                foreach (Atom a in m.Atoms) U += a.U;
                H += U;
            }
            return H;
        }
        public string Solve()
        {
            Console.WriteLine(0);
            Data.WriteDate();
            foreach (Molecule m in Molecules)
            {
                m.CoofTenzor();
                m.KineticMoment();
            }
            for (int t = 0; t < T; t++)
            {
                foreach (Molecule m in Molecules)
                {
                    foreach (Atom a in m.Atoms)
                    {
                        a.x0 = a.x;
                        a.y0 = a.y;
                        a.z0 = a.z;
                        a.xA = a.x + m.x;
                        a.yA = a.y + m.y;
                        a.zA = a.z + m.z;
                    }
                    m.x0 = m.x;
                    m.y0 = m.y;
                    m.z0 = m.z;
                    m.u0 = m.u;
                    m.v0 = m.v;
                    m.w0 = m.w;
                    m.Kx0 = m.Kx;
                    m.Ky0 = m.Ky;
                    m.Kz0 = m.Kz;
                }
                //////1
                UVW();
                foreach (Molecule m in Molecules)
                {
                    foreach (Atom a in m.Atoms)
                    {
                        a.VelocityStep();
                        a.xb1 = a.u;
                        a.yb1 = a.v;
                        a.zb1 = a.w;
                    }
                    m.MomentForces();
                    m.ub1 = m.ub;
                    m.vb1 = m.vb;
                    m.wb1 = m.wb;
                    m.xb1 = m.u;
                    m.yb1 = m.v;
                    m.zb1 = m.w;
                    m.Lx1 = m.Lx;
                    m.Ly1 = m.Ly;
                    m.Lz1 = m.Lz;
                    //
                    m.x = m.x0 + dt / 2.0 * m.xb1;
                    m.y = m.y0 + dt / 2.0 * m.yb1;
                    m.z = m.z0 + dt / 2.0 * m.zb1;
                    m.u = m.u0 + dt / 2.0 * m.ub1;
                    m.v = m.v0 + dt / 2.0 * m.vb1;
                    m.w = m.w0 + dt / 2.0 * m.wb1;
                    foreach (Atom a in m.Atoms)
                    {
                        a.x = a.x0 + dt / 2.0 * a.xb1;
                        a.y = a.y0 + dt / 2.0 * a.yb1;
                        a.z = a.z0 + dt / 2.0 * a.zb1;
                        a.xA = a.x + m.x;
                        a.yA = a.y + m.y;
                        a.zA = a.z + m.z;
                    }
                    m.Kx = m.Kx0 + dt / 2.0 * m.Lx1;
                    m.Ky = m.Ky0 + dt / 2.0 * m.Ly1;
                    m.Kz = m.Kz0 + dt / 2.0 * m.Lz1;
                    m.CoofTenzor();
                    m.AngleVelocity();
                }
                //////2
                UVW();
                foreach (Molecule m in Molecules)
                {
                    foreach (Atom a in m.Atoms)
                    {
                        a.VelocityStep();
                        a.xb2 = a.u;
                        a.yb2 = a.v;
                        a.zb2 = a.w;
                    }
                    m.MomentForces();
                    m.ub2 = m.ub;
                    m.vb2 = m.vb;
                    m.wb2 = m.wb;
                    m.xb2 = m.u;
                    m.yb2 = m.v;
                    m.zb2 = m.w;
                    m.Lx2 = m.Lx;
                    m.Ly2 = m.Ly;
                    m.Lz2 = m.Lz;
                    //
                    m.x = m.x0 + dt / 2.0 * m.xb2;
                    m.y = m.y0 + dt / 2.0 * m.yb2;
                    m.z = m.z0 + dt / 2.0 * m.zb2;
                    m.u = m.u0 + dt / 2.0 * m.ub2;
                    m.v = m.v0 + dt / 2.0 * m.vb2;
                    m.w = m.w0 + dt / 2.0 * m.wb2;
                    foreach (Atom a in m.Atoms)
                    {
                        a.x = a.x0 + dt / 2.0 * a.xb2;
                        a.y = a.y0 + dt / 2.0 * a.yb2;
                        a.z = a.z0 + dt / 2.0 * a.zb2;
                        a.xA = a.x + m.x;
                        a.yA = a.y + m.y;
                        a.zA = a.z + m.z;
                    }
                    m.Kx = m.Kx0 + dt / 2.0 * m.Lx2;
                    m.Ky = m.Ky0 + dt / 2.0 * m.Ly2;
                    m.Kz = m.Kz0 + dt / 2.0 * m.Lz2;
                    m.CoofTenzor();
                    m.AngleVelocity();
                }
                //////3
                UVW();
                foreach (Molecule m in Molecules)
                {
                    foreach (Atom a in m.Atoms)
                    {
                        a.VelocityStep();
                        a.xb3 = a.u;
                        a.yb3 = a.v;
                        a.zb3 = a.w;
                    }
                    m.MomentForces();
                    m.ub3 = m.ub;
                    m.vb3 = m.vb;
                    m.wb3 = m.wb;
                    m.xb3 = m.u;
                    m.yb3 = m.v;
                    m.zb3 = m.w;
                    m.Lx3 = m.Lx;
                    m.Ly3 = m.Ly;
                    m.Lz3 = m.Lz;
                    //
                    m.x = m.x0 + dt * m.xb3;
                    m.y = m.y0 + dt * m.yb3;
                    m.z = m.z0 + dt * m.zb3;
                    m.u = m.u0 + dt * m.ub3;
                    m.v = m.v0 + dt * m.vb3;
                    m.w = m.w0 + dt * m.wb3;
                    foreach (Atom a in m.Atoms)
                    {
                        a.x = a.x0 + dt * a.xb3;
                        a.y = a.y0 + dt * a.yb3;
                        a.z = a.z0 + dt * a.zb3;
                        a.xA = a.x + m.x;
                        a.yA = a.y + m.y;
                        a.zA = a.z + m.z;
                    }
                    m.Kx = m.Kx0 + dt * m.Lx3;
                    m.Ky = m.Ky0 + dt * m.Ly3;
                    m.Kz = m.Kz0 + dt * m.Lz3;
                    m.CoofTenzor();
                    m.AngleVelocity();
                }
                //////4
                UVW();
                foreach (Molecule m in Molecules)
                {
                    foreach (Atom a in m.Atoms)
                    {
                        a.VelocityStep();
                        a.xb4 = a.u;
                        a.yb4 = a.v;
                        a.zb4 = a.w;
                    }
                    m.MomentForces();
                    m.ub4 = m.ub;
                    m.vb4 = m.vb;
                    m.wb4 = m.wb;
                    m.xb4 = m.u;
                    m.yb4 = m.v;
                    m.zb4 = m.w;
                    m.Lx4 = m.Lx;
                    m.Ly4 = m.Ly;
                    m.Lz4 = m.Lz;
                    //
                    m.x = m.x0 + dt / 6.0 * (m.xb1 + 2.0 * m.xb2 + 2.0 * m.xb3 + m.xb4);
                    m.y = m.y0 + dt / 6.0 * (m.yb1 + 2.0 * m.yb2 + 2.0 * m.yb3 + m.yb4);
                    m.z = m.z0 + dt / 6.0 * (m.zb1 + 2.0 * m.zb2 + 2.0 * m.zb3 + m.zb4);
                    m.u = m.u0 + dt / 6.0 * (m.ub1 + 2.0 * m.ub2 + 2.0 * m.ub3 + m.ub4);
                    m.v = m.v0 + dt / 6.0 * (m.vb1 + 2.0 * m.vb2 + 2.0 * m.vb3 + m.vb4);
                    m.w = m.w0 + dt / 6.0 * (m.wb1 + 2.0 * m.wb2 + 2.0 * m.wb3 + m.wb4);
                    foreach (Atom a in m.Atoms)
                    {
                        a.x = a.x0 + dt / 6.0 * (a.xb1 + 2.0 * a.xb2 + 2.0 * a.xb3 + a.xb4);
                        a.y = a.y0 + dt / 6.0 * (a.yb1 + 2.0 * a.yb2 + 2.0 * a.yb3 + a.yb4);
                        a.z = a.z0 + dt / 6.0 * (a.zb1 + 2.0 * a.zb2 + 2.0 * a.zb3 + a.zb4);
                    }
                    m.Kx = m.Kx0 + dt / 6.0 * (m.Lx1 + 2.0 * m.Lx2 + 2.0 * m.Lx3 + m.Lx4);
                    m.Ky = m.Ky0 + dt / 6.0 * (m.Ly1 + 2.0 * m.Ly2 + 2.0 * m.Ly3 + m.Ly4);
                    m.Kz = m.Kz0 + dt / 6.0 * (m.Lz1 + 2.0 * m.Lz2 + 2.0 * m.Lz3 + m.Lz4);
                    m.CoofTenzor();
                    m.AngleVelocity();
                }
                Data.WriteDate();
                if (T >= 1000)
                {
                    if (t % (T / 1000) == 0)
                    {
                        Data.WriteEnergy(Balance());
                        Console.SetCursorPosition(0, Console.CursorTop - 1);
                        Console.Write("\r" + new string(' ', Console.BufferWidth) + "\r");
                        Console.WriteLine(((t * 1000.0 / T + 1)/10.0).ToString()+"%");
                    }
                }
                else Console.WriteLine(t.ToString());
            }
            return "Completed!!!";
        }
    }
    static class Data
    {
        static int T; static double dt;
        public static void SetTdt(int T, double dt)
        {
            Data.T = T; Data.dt = dt;
        }
        private static readonly List<Molecule> Molecules = new List<Molecule>();
        private static readonly List<StreamWriter> Lsw = new List<StreamWriter>();
        static StreamWriter sw1;
        public static void ReadFiles(string s, Molecule M)
        {
            Molecules.Add(M);
            string RootDirectory = Directory.GetCurrentDirectory();
            s = RootDirectory + "\\Data\\" + s;
            string[] Data = File.ReadAllLines(s);
            for (int i = 0; i < Data.Length; i++)
            {
                if (Data[i] == "Carbon")
                {
                    i++;
                    int N = int.Parse(Data[i]);
                    i++;
                    for (int j = 0; j < N; j++)
                    {
                        string[] xyz = Data[i + j].Split(new char[] { ' ' });
                        M.AddAtom(new Carbon(double.Parse(xyz[0], CultureInfo.InvariantCulture), double.Parse(xyz[1], CultureInfo.InvariantCulture), double.Parse(xyz[2], CultureInfo.InvariantCulture), M));
                    }
                    i++;
                }
            }
        }
        public static void OpenFiles(string str)
        {
            string RootDirectory = Directory.GetCurrentDirectory();
            str = RootDirectory + "\\Results\\" + str;
            Directory.CreateDirectory(str);
            FileStream fs0 = new FileStream(str + "\\Description.txt", FileMode.Create, FileAccess.Write);
            StreamWriter sw0 = new StreamWriter(fs0);
            FileStream fs1 = new FileStream(str + "\\Energy.txt", FileMode.Create, FileAccess.Write);
            sw1 = new StreamWriter(fs1);
            sw0.WriteLine("Description:");
            sw0.WriteLine("T = " + Data.T.ToString() + ", dt = " + Data.dt.ToString() + ";");
            for (int i = 0; i < Molecules.Count; i++)
            {
                sw0.WriteLine(Molecules[i].name);
                sw0.WriteLine("   Initiate Data:");
                sw0.WriteLine("   x = " + Molecules[i].xM0.ToString() + ", y = " + Molecules[i].yM0.ToString() + ", z = " + Molecules[i].zM0.ToString() + ";");
                sw0.WriteLine("   u = " + Molecules[i].uM0.ToString() + ", v = " + Molecules[i].vM0.ToString() + ", w = " + Molecules[i].wM0.ToString() + ";");
                sw0.WriteLine("   OmegaX = " + Molecules[i].OmegaX0.ToString() + ", OmegaY = " + Molecules[i].OmegaY0.ToString() + ", OmegaZ = " + Molecules[i].OmegaZ0.ToString() + ";");
                FileStream fs = new FileStream(str + "\\" + Molecules[i].name + ".txt", FileMode.Create, FileAccess.Write);
                Lsw.Add(new StreamWriter(fs));
            }
            sw0.Close();
        }
        public static void WriteDate()
        {
            string s;
            for (int i = 0; i < Molecules.Count; i++)
            {
                s = "";
                s += Molecules[i].x.ToString() + ' ' + Molecules[i].y.ToString() + ' ' + Molecules[i].z.ToString() + ' ';
                s += Molecules[i].u.ToString() + ' ' + Molecules[i].v.ToString() + ' ' + Molecules[i].w.ToString();
                Lsw[i].WriteLine(s.Replace(',', '.'));
            }
        }
        public static void WriteEnergy(double Energy)
        {
            sw1.WriteLine(Energy.ToString().Replace(',', '.'));
        }
        public static void CloseFiles()
        {
            foreach (StreamWriter sw in Lsw)
            {
                sw.Close();
            }
            sw1.Close();
            Console.WriteLine("Files writing completed");
        }
    }
    static class Constant
    {
        public static double KB = 1.380649 * Math.Pow(10.0, -23);
        public static double atomMass = 1.6605390666 * Math.Pow(10.0, -27);
    }
    internal static class Program
    {
        static void Main()
        {
            //Initiate
            Molecule Molecule1 = new Molecule("Первая углеродная нанотрубка");
            Molecule Molecule2 = new Molecule("Вторая углеродная нанотрубка");
            const string Nanotube = "NanoTube.txt";
            const string Fullerene = "Fullerene.txt";
            Data.ReadFiles(Nanotube, Molecule1);
            Data.ReadFiles(Nanotube, Molecule2);
            Molecule1.CMToZero();
            Molecule2.CMToZero();
            Molecule1.InitiateSpaceCoor(-2.5, -1.0, 0.0); Molecule1.InitiateVelocityCoor(12.0, 0.0, 0.0); Molecule1.InitiateAngleVelocityCoor(0.0, 0.0, 1299.7);
            Molecule2.InitiateSpaceCoor(2.5, 1.0, 0.0); Molecule2.InitiateVelocityCoor(-12.0, 0.0, 0.0); Molecule2.InitiateAngleVelocityCoor(0.0, 0.0, 1299.7);
            Solver Solver = new Solver(60000, Math.Pow(10.0, -5.0));
            Solver.AddMolecule(Molecule1);
            Solver.AddMolecule(Molecule2);
            //Data
            Data.OpenFiles("ResultsNanotube26.06.24");
            //Solver
            Console.WriteLine(Solver.Solve());
            Data.CloseFiles();
            //EndProgramm
            Console.ReadKey();
        }
    }
}
