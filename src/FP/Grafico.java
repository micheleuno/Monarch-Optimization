package FP;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import MCDP.model.MCDPData;


@SuppressWarnings("serial")
public class Grafico extends javax.swing.JFrame {
    
    JFreeChart Graphic;
    XYSeries series = new XYSeries("Resultado");
    XYDataset Datos = new XYSeriesCollection(series);
    
    private int nFlores;
    private int nPartes;
    private int nCeldas;
    private int nIteraciones;
    private long mejorFitness;
    private int maxMaquinas;
    private String nProblema;
    
    public Grafico(int vector[], MCDPData data, int numberIteration, int numberPoblation) {
        initComponents();
        nFlores = numberPoblation;
        nIteraciones = numberIteration;
        nProblema = data.getIdentificator();
        nPartes = data.getP();
        nCeldas = data.getC();
        maxMaquinas = data.getMmax();
        mejorFitness = data.getBestSGlobal();

        int mejorFitness = 10000;
        for(int i = 0; i < vector.length; i++){
            if(vector[i] < mejorFitness) mejorFitness = vector[i];
            series.add((i+1), vector[i]);
        }
        Graphic = ChartFactory.createXYLineChart("Problema "+nProblema, "Ejecuci�n", "Fitness", Datos, PlotOrientation.VERTICAL, false, true, true);
 
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        boton_mostrar = new javax.swing.JButton();
        jLabel1 = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        jLabel3 = new javax.swing.JLabel();
        jLabel4 = new javax.swing.JLabel();
        jLabel5 = new javax.swing.JLabel();
        jLabel6 = new javax.swing.JLabel();
        jLabel9 = new javax.swing.JLabel();
        res_problema = new javax.swing.JLabel();
        res_celdas = new javax.swing.JLabel();
        res_partes = new javax.swing.JLabel();
        res_max = new javax.swing.JLabel();
        res_it = new javax.swing.JLabel();
        res_bats = new javax.swing.JLabel();
        res_best = new javax.swing.JLabel();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);

        boton_mostrar.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        boton_mostrar.setText("Graficar");
        boton_mostrar.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                boton_mostrarActionPerformed(evt);
            }
        });

        jLabel1.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        jLabel1.setText("N� Celdas:");

        jLabel2.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        jLabel2.setText("Problema N�:");

        jLabel3.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        jLabel3.setHorizontalAlignment(javax.swing.SwingConstants.LEFT);
        jLabel3.setText("N� Partes:");

        jLabel4.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        jLabel4.setText("Max M�quinas:");

        jLabel5.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        jLabel5.setText("Iteraciones:");

        jLabel6.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        jLabel6.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
        jLabel6.setText("N� Flores:");

        jLabel9.setFont(new java.awt.Font("Tahoma", 0, 24)); // NOI18N
        jLabel9.setForeground(new java.awt.Color(0, 153, 0));
        jLabel9.setText("Mejor Fitness: ");

        res_problema.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N

        res_celdas.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N

        res_partes.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N

        res_max.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N

        res_it.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N

        res_bats.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N

        res_best.setFont(new java.awt.Font("Tahoma", 1, 24)); // NOI18N
        res_best.setForeground(new java.awt.Color(0, 153, 0));

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(19, 19, 19)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(jLabel3)
                    .addComponent(jLabel1)
                    .addComponent(jLabel4)
                    .addComponent(jLabel2))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(res_problema)
                            .addComponent(res_celdas))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 78, Short.MAX_VALUE)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(jLabel6, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jLabel5, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(res_it))
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                                .addGap(10, 10, 10)
                                .addComponent(res_bats)))
                        .addGap(75, 75, 75))
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(res_partes)
                            .addComponent(res_max))
                        .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(136, 136, 136)
                        .addComponent(boton_mostrar))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(68, 68, 68)
                        .addComponent(jLabel9)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(res_best)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addGap(19, 19, 19)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel2)
                    .addComponent(jLabel5)
                    .addComponent(res_problema)
                    .addComponent(res_it))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(jLabel6)
                    .addComponent(res_celdas)
                    .addComponent(res_bats))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel3)
                    .addComponent(res_partes))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel4)
                    .addComponent(res_max))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 28, Short.MAX_VALUE)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel9)
                    .addComponent(res_best))
                .addGap(26, 26, 26)
                .addComponent(boton_mostrar)
                .addGap(53, 53, 53))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void boton_mostrarActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_boton_mostrarActionPerformed
        ChartPanel Panel = new ChartPanel(Graphic);
        JFrame Ventana = new JFrame("Resultados");
        Ventana.getContentPane().add(Panel);
        Ventana.pack();
        Ventana.setVisible(true);
        Ventana.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        this.res_bats.setText(Integer.toString(nFlores));
        this.res_it.setText(Integer.toString(nIteraciones));
        this.res_celdas.setText(Integer.toString(nCeldas));
        this.res_partes.setText(Integer.toString(nPartes));
        this.res_max.setText(Integer.toString(maxMaquinas));
        this.res_best.setText(Long.toString(mejorFitness));
        this.res_problema.setText(nProblema);
    }//GEN-LAST:event_boton_mostrarActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton boton_mostrar;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel9;
    private javax.swing.JLabel res_bats;
    private javax.swing.JLabel res_best;
    private javax.swing.JLabel res_celdas;
    private javax.swing.JLabel res_it;
    private javax.swing.JLabel res_max;
    private javax.swing.JLabel res_partes;
    private javax.swing.JLabel res_problema;
    // End of variables declaration//GEN-END:variables
}
