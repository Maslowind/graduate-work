using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace Graduate_work
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        public class Result
        {
            public double R { get; set; }
            public double S_tt { get; set; }
            public double S_st { get; set; }
            public double S_ts { get; set; }
            public double S_ss { get; set; }
            public double sigma_rr { get; set; }
            public double u_r { get; set; }
            public double sigma_fifi { get; set; }
        }
        double nu, mu, alpha, beta, gamma;
        double R0 = 1, R1, R, h;
        double  tau_coef_R1, A_coef_R1;
        List<Result> staticResults = new List<Result>();
        List<Result> results = new List<Result>();
        List<Result> results_lin = new List<Result>();
        List<Result> results_quad = new List<Result>();
        List<Result> results_cube = new List<Result>();
        Func<double, double> G_proiz;
        List<Func<double, double>> G_proiz_list = new List<Func<double, double>>();
        double S_tt_static(double R)
        {
            return Math.Pow(R / R0, alpha) * (mu * (1 - gamma / alpha) - (mu + alpha) * Math.Pow(R / R0, beta)) / beta;
        }
        double S_st_static(double R)
        {
            return mu * Math.Pow(R / R0, alpha) * (Math.Pow(R / R0, beta) - 1) / beta;
        }
        double S_ts_static(double R)
        {
            return Math.Pow(R / R0, alpha) * (Math.Pow(R / R0, beta) - 1) / (beta * (1 - nu));
        }
        double S_ss_static(double R)
        {
            return Math.Pow(R / R0, alpha) * (-1 * (mu + alpha) + mu * (1 - gamma / alpha) * Math.Pow(R / R0, beta)) / beta;
        }

        double G_static(double R)
        {
            return Math.Pow(R, gamma);
        }
        double G_lin(double R)
        {
            return 2 * R + 20;
        }
        double G_quad(double R)
        {
            return 10 * R * R + 5 * R + 43;
        }
        double G_cube(double R)
        {
            return R*R*R + 20*R*R + 7;
        }
        



        double ln_G_static_proiz(double R)
        {
            return gamma / R;
        }
        double ln_G_lin_proiz(double R)
        {
            return 2.0 / G_lin(R);
        }
        double ln_G_quad_proiz(double R)
        {
            return 1.0 * (20 * R + 5) / G_quad(R);
        }
        double ln_G_cube_proiz(double R)
        {
            return 1.0*(3 * R*R+40*R) / G_cube(R);
        }


        double sigma_rr_func(double R, List<Result> list)
        {
            int id = list.FindIndex(x => Math.Round(x.R, 3) == Math.Round(R, 3));
            return (list[id].S_tt - A_coef_R1 * list[id].S_ts) * tau_coef_R1;
        }
        double u_r_func(double R, List<Result> list)
        {
            int id = list.FindIndex(x => Math.Round(x.R, 3) == Math.Round(R, 3));
            return (list[id].S_st - A_coef_R1 * list[id].S_ss) * tau_coef_R1;
        }
        double sigma_fifi_func(double R, List<Result> list)
        {
            int id = list.FindIndex(x => Math.Round(x.R, 3) == Math.Round(R, 3));
            return (nu * list[id].sigma_rr+ list[id].u_r)/(1+nu);
        }




        double sigma_rr_proiz(double R, double sigma_rr, double s)
        {
            double result = s / ((1 - nu) * R) - (1 - 2 * nu) * sigma_rr / ((1 - nu) * R);
            return result;
        }
        double s_proiz(double R, double s, double sigma_rr)
        {
            double result = (1 - 2 * nu) * sigma_rr / ((1 - nu) * R) - s / ((1 - nu) * R) + s * G_proiz(R);
            return result;
        }

        void  getCoef(List<Result> list)
        {
            int id = list.FindIndex(x => Math.Round(x.R,2) == Math.Round(R1,2));
            double a = list[id].S_tt;
            double b = list[id].S_ts;
            double c = list[id].S_st;
            double d = list[id].S_ss;
            tau_coef_R1 = 1;
            A_coef_R1 = 1.0*c/d;
        }

        private double methodRungeKutta(double x, double y1, double y2, double h, bool leftOrRight, Func<double, double, double, double> func)
        {
            double k1, k2, k3, k4, koef;
            if (leftOrRight) koef = 1;
            else koef = -1;
            h = h * koef;
            k1 = h * func(x, y1, y2);
            k2 = h * func(x + h / 2, y1 + koef * k1 / 2, y2);
            k3 = h * func(x + h / 2, y1 + koef * k2 / 2, y2);
            k4 = h * func(x + h, y1 + koef * k3, y2);
            double result = y1 + koef * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            return result;
        }

        void addFuncToList()
        {
            G_proiz_list.Add(new Func<double, double>(ln_G_static_proiz));
            G_proiz_list.Add(new Func<double, double>(ln_G_lin_proiz));
            G_proiz_list.Add(new Func<double, double>(ln_G_quad_proiz));
            G_proiz_list.Add(new Func<double, double>(ln_G_cube_proiz));
        }

        void addResults(int i, double R, double sigma_rr_1, double s_0, double sigma_rr_0, double s_1,ref List<Result> list)
        {

            list.Add(new Result());            
            list[i].R = R;
            list[i].S_tt = sigma_rr_1;
            list[i].S_st = s_0;
            list[i].S_ts = sigma_rr_0;
            list[i].S_ss = s_1;
        }
        void addFinalResults( ref List<Result> list)
        {
            for (int j = 0; j < list.Count; j++)
            {
                list[j].sigma_rr = sigma_rr_func(list[j].R, list);
                list[j].u_r = u_r_func(list[j].R, list);
                list[j].sigma_fifi = sigma_fifi_func(list[j].R, list);
            }
        }

        void calculateEverything(ref List<Result> list, string name, Color color)
        {
            R = R0;
            int i = 1;
            h = 0.05;
            double sigma_rr_1 = 1, s_0 = 0, sigma_rr_0 = 0, s_1 = 1;
            addResults(0, R, sigma_rr_1, s_0, sigma_rr_0, s_1, ref list);
            
            while (R <= R1)
            {
                R += h;
                double sigma_rr_1_2 = methodRungeKutta(R, sigma_rr_1, s_0, h, true, sigma_rr_proiz);
                double s_0_2 = methodRungeKutta(R, s_0, sigma_rr_1, h, true, s_proiz);
                double sigma_rr_0_2 = methodRungeKutta(R, sigma_rr_0, s_1, h, true, sigma_rr_proiz);
                double s_1_2 = methodRungeKutta(R, s_1, sigma_rr_0, h, true, s_proiz);
                sigma_rr_1 = sigma_rr_1_2; s_0 = s_0_2; sigma_rr_0 = sigma_rr_0_2; s_1 = s_1_2;

                addResults(i, R, sigma_rr_1, s_0, sigma_rr_0, s_1, ref list);
                i++;
            }
            getCoef(list);
            addFinalResults(ref list);
            addToChart(list, name, color);
        }
        void addToChart(List<Result> list, string name, Color color)
        {
            chart1.ChartAreas[0].AxisX.Maximum = R1; chart1.ChartAreas[0].AxisX.Minimum = R0;
            chart2.ChartAreas[0].AxisX.Maximum = R1; chart2.ChartAreas[0].AxisX.Minimum = R0;
            chart3.ChartAreas[0].AxisX.Maximum = R1; chart3.ChartAreas[0].AxisX.Minimum = R0;
            Random rnd = new Random();
            Series func_sigma_rr = new Series(name);
            Series func_u_r = new Series(name);
            Series func_sigma_fifi = new Series(name);
            func_sigma_rr.ChartType = SeriesChartType.Line;
            func_sigma_rr.Color = color;
            func_u_r.ChartType = SeriesChartType.Line;
            func_u_r.Color = color;
            func_sigma_fifi.ChartType = SeriesChartType.Line;
            func_sigma_fifi.Color = color;
            for (int i = 0; i < list.Count; i++)
            {
                func_sigma_rr.Points.AddXY(list[i].R, list[i].sigma_rr);
                func_u_r.Points.AddXY(list[i].R, list[i].u_r);
                func_sigma_fifi.Points.AddXY(list[i].R, list[i].sigma_fifi);
            }
            chart1.Series.Add(func_sigma_rr);
            chart2.Series.Add(func_u_r);
            chart3.Series.Add(func_sigma_fifi);
        }
        void Reset ()
        {
            staticResults.Clear();
            results.Clear();
            results_lin.Clear();
            results_quad.Clear();
            results_cube.Clear();
            listBox1.Items.Clear();
            listBox1.Items.Add("R\t\nS_ττ\t\nS_sτ\t\nS_τs\t\nS_ss");
            listBox2.Items.Clear();
            listBox2.Items.Add("R\t\nS_ττ\t\nS_sτ\t\nS_τs\t\nS_ss");
            while (chart1.Series.Count > 0) { chart1.Series.RemoveAt(0); }
            while (chart2.Series.Count > 0) { chart2.Series.RemoveAt(0); }
            while (chart3.Series.Count > 0) { chart3.Series.RemoveAt(0); }
        }

        private void button1_Click(object sender, EventArgs e)
        {
            Reset();            
            R1 = Convert.ToDouble(textBox2.Text);
            nu = Convert.ToDouble(textBox3.Text);
            gamma = Convert.ToDouble(textBox1.Text);
            mu = 1.0 * (1 - 2 * nu) / (1 - nu);
            beta = Math.Sqrt(gamma * gamma - (4.0 * nu * gamma) / (1 - nu) + 4);
            alpha = 0.5 * (gamma - beta) - 1;

            addFuncToList();
            G_proiz = G_proiz_list[0];
            calculateEverything(ref results, "Степенева", Color.FromName("Red"));
            G_proiz = G_proiz_list[1];
            calculateEverything(ref results_lin, "Лінійна", Color.FromName("Blue"));
            G_proiz = G_proiz_list[2];
            calculateEverything(ref results_quad, "Квадратична", Color.FromName("Green"));
            G_proiz = G_proiz_list[3];
            calculateEverything(ref results_cube, "Кубічна", Color.FromName("Black"));
            R = R0;
            h = 0.05;
            for (int i = 0; i < results.Count; i++)
            {

                addResults(i, R, S_tt_static(R), S_st_static(R), S_ts_static(R), S_ss_static(R), ref staticResults);
                listBox1.Items.Add(Math.Round(results[i].R, 3) + "\t\n" + Math.Round(results[i].S_tt, 3) + "\t\n" + Math.Round(results[i].S_st, 3) + "\t\n" + Math.Round(results[i].S_ts, 3) + "\t\n" + Math.Round(results[i].S_ss, 3) + "\t\n");
                listBox2.Items.Add(Math.Round(staticResults[i].R, 3) + "\t\n" + Math.Round(staticResults[i].S_tt, 3) + "\t\n" + Math.Round(staticResults[i].S_st, 3) + "\t\n" + Math.Round(staticResults[i].S_ts, 3) + "\t\n" + Math.Round(staticResults[i].S_ss, 3) + "\t\n");
                
                R += h;
            }
            
        }
    }
}
