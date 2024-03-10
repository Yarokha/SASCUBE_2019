#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    gridLayout = new QGridLayout(ui->widgetPlot);
    gridLayout2D = new QGridLayout(ui->widget2DPlot);

    ui->NormalBox->addItem("X");
    ui->NormalBox->addItem("Y");
    ui->NormalBox->addItem("Z");

}

MainWindow::~MainWindow()
{
    delete ui;
}


double MainWindow::mult(std::vector<double>& data1, std::vector<double>& data2) {
    double sum = 0;
    for (int i = 0; i < 3; ++i)
        sum += data1[i] * data2[i];
    return sum;
}

void MainWindow::Rotate(std::vector<double>& data, double a) {
    std::vector<std::vector<double>> rot_mat = { { 1,0,0 },{ 0,cos(a),-sin(a) },{ 0,sin(a),cos(a) } };
    std::vector<double> temp(0);
    for (int i = 0; i < 3; ++i)
        temp.push_back(mult(rot_mat[i], data));

    for (int i = 0; i < 3; ++i)
        data[i] = temp[i];
}

void MainWindow::Output_File_Write_(QTextEdit *log){
    QFile file(OutputFilePath_);
    if (file.open(QFile::WriteOnly|QFile::Truncate) )
    {
           QTextStream cout(&file);
           cout<<log->toPlainText();
    }
    file.close();
}

void MainWindow::Log_state_(QString Process, QTextEdit *log, QProgressBar *Bar, int percent_state){
    std::time_t time = std::time(nullptr);
    if (percent_state == 0)
        log->clear();

    log->append(Process);
    log->append(std::asctime(std::localtime(&time)));
    Bar->setValue(percent_state);
    Output_File_Write_(log);
}

void MainWindow::on_PDBchooseButon_clicked()
{
    ui->ProtInfo->clear();
    QString Path = QFileDialog::getOpenFileName(this,
        tr("Open Protein Data Bank File"), "/", tr("PDB Files (*.pdb)"));
    ui->PDBfilePath->setText(Path);

    QFileInfo File(Path);

    ProjectPath_ = Path.left(Path.length()-4);

    if (!QDir(ProjectPath_).exists()){
            QDir().mkdir(ProjectPath_);
    }

    PipedFilePath_ = ProjectPath_ + "/piped.txt";
    QFileInfo check_pipedfile(PipedFilePath_);
    if (check_pipedfile.exists())
        ui->PipedFilePath->setText(PipedFilePath_);
    else
        ui->PipedFilePath->setText("");

    QString CurveFilePath_ = ProjectPath_ + "/curve.txt";
    QFileInfo check_curvefile(CurveFilePath_);
    if (check_curvefile.exists())
        ui->PlotButton->setEnabled(true);
    else
        ui->PlotButton->setEnabled(false);

    QString ProteinFilePath_ = ProjectPath_ + "/Molecular.txt";
    QFileInfo check_proteinfile(ProteinFilePath_);
    if (check_proteinfile.exists()){
        ui->Plot2DButton->setEnabled(true);
        ui->ButtonDown->setEnabled(true);
        ui->ButtonUp->setEnabled(true);
    }
    else{
        ui->Plot2DButton->setEnabled(false);
        ui->ButtonDown->setEnabled(false);
        ui->ButtonUp->setEnabled(false);
    }


    QFile Param(ProjectPath_+"/ProjectParameters.txt");
    if (Param.open(QFile::ReadOnly|QFile::Truncate)){
        ui->ProtInfo->clear();
        QTextStream stream(&Param);
        while (!stream.atEnd()){
            ui->ProtInfo->append(stream.readLine());
        }
    }
    Param.close();

    if(Path!="")
        ui->StartButton->setEnabled(true);
    else{
        ui->StartButton->setEnabled(false);
        ui->PlotButton->setEnabled(false);
    }
    OutputFilePath_ = ProjectPath_+"/output.txt";
    ui->statusBar->showMessage("", 0);
    ui->LogEdit->clear();
}

void MainWindow::on_PipedCheckBox_clicked()
{
    if (ui->PipedCheckBox->checkState()==false){
        ui->PipedFilePath->setEnabled(false);
        ui->PipedChooseButon->setEnabled(false);

        ui->rho0->setEnabled(true);
        ui->CubeRib->setEnabled(true);
        ui->SolventRadius->setEnabled(true);

        ui->label_9->setEnabled(true);
        ui->label_7->setEnabled(true);
        ui->label_2->setEnabled(true);
    }
    else{
        ui->PipedFilePath->setEnabled(true);
        ui->PipedChooseButon->setEnabled(true);


        ui->rho0->setEnabled(false);
        ui->CubeRib->setEnabled(false);
        ui->SolventRadius->setEnabled(false);

        ui->label_9->setEnabled(false);
        ui->label_7->setEnabled(false);
        ui->label_2->setEnabled(false);
    }
}

void MainWindow::on_PipedChooseButon_clicked()
{
    PipedFilePath_ = QFileDialog::getOpenFileName(this,
        tr("Open File of Parallelepiped"), "/", tr("piped.txt"));
    ui->PipedFilePath->setText(PipedFilePath_);
}

void MainWindow::on_actionExit_triggered()
{
    exit(0);
}

void MainWindow::on_actionSettings_triggered()
{
    Dialog = new PDBSettings(this);
    Dialog->show();
}

void MainWindow::on_ExitButton_clicked()
{
    exit(0);
}

void MainWindow::on_StartButton_clicked()
{
//    Names
    QString
    I_p_str = "Input parameters",
    c_e_str = "Cube edge = ",
    S_R_str = "Solvent Radius = ",
    S_E_D_str = "Solvent Electronic Density = ",
    Q_m_str = "Qmin = ",
    Q_M_str = "Qmax = ",
    N_Q_str = "Number of Q = ",
    N_p_str = "Nphi = ",
    N_z_str = "Nz = ",
    C_p_str = "Calculated properties:",
    H_v_str = "Hydrated volume of the protein = ",
    N_C_h_s_str = "Number of cubes of the hydrated surface of the protein = ",
    M_v_str = "Molecular volume of the protein = ",
    N_c_m_s_str = "Number of cubes of the molecular surface of the protein = ",
    C_str = "Chi = ",
    G_R_str = "Gyration Radius = ";
//-----------------------------------------------------------------------------------------------


    ui->statusBar->showMessage("Calculation in progress", 0);
    std::time_t Start_TIME = std::time(nullptr);
    ui->PlotButton->setEnabled(false);
    ui->StartButton->setEnabled(false);
    ui->Plot2DButton->setEnabled(false);
    ui->ButtonDown->setEnabled(false);
    ui->ButtonUp->setEnabled(false);
    Log_state_("Start", ui->LogEdit, ui->progressBar, 0);
    double solvent = (ui->SolventRadius->text()).toDouble(),
            a = (ui->CubeRib->text()).toDouble(),
            //a = (ui->CubeBox->currentText()).toDouble(),
            Q_max = (ui->Qmax->text()).toDouble(),
            Q_min = (ui->Qmin->text()).toDouble(),
            rho_0 = (ui->rho0->text()).toDouble(),
            delta = 0.015;// (ui->Delta->text()).toDouble();
    int Q_n = (ui->Qn->text()).toInt();

    freopen((ProjectPath_.toStdString() + "/log.txt").c_str(), "w", stdout);

    std::vector< std::string > sett;
    std::string temp;

    Log_state_("Reading PDB File", ui->LogEdit, ui->progressBar, 5);
    std::ifstream set("Settings.txt" );
    while(std::getline(set, temp)){
        sett.push_back(temp);
    }



    SAXSParser ProteinTab(ui->PDBfilePath->text().toStdString(), sett[1], sett[0]);
    ProteinTab.SaveDataToFile(true, true, (ui->checkBoxHET->checkState()), ui->checkBoxHOH->checkState(), (ui->checkBoxHEM->checkState()), ProjectPath_.toStdString());
    std::vector < std::vector<double> > PDB_table = ProteinTab.AtomTable(false, true, (ui->checkBoxHET->checkState()), ui->checkBoxHOH->checkState(), (ui->checkBoxHEM->checkState()));
    SAScore Protein(PDB_table, solvent, a, delta);
    Protein.big_piped = ProteinTab.big_piped;
    Protein.ProjectFolder(ProjectPath_.toStdString());
    Protein.SetQSpace(ui->NzLabel->text().toInt(), ui->NphiLabel->text().toInt());
//    if (!ui->PipedCheckBox->checkState()){

      if (ui->GeomCheckBox->checkState()){
        Log_state_("Building Geometry", ui->LogEdit, ui->progressBar, 30);
        Protein.CoreStart();
      }
//    }
//    else{
//        Log_state_("Loading Parallelepipeds", ui->LogEdit, ui->progressBar, 30);
//        Protein.PipedLoad((ui->PipedFilePath->text()).toStdString());
//    }

    Log_state_("Calculating Intensity", ui->LogEdit, ui->progressBar, 50);
    Protein.CalcIntensityBigPiped(Q_min, Q_max, Q_n, rho_0);

//    if (!ui->ExperCheckBox->checkState())
//        Protein.CalcIntensity(Q_min, Q_max, Q_n, rho_0);
//    else{
//        Protein.CalculatingIntChi(Ie_, rho_0);
//    }

    Log_state_("Finish", ui->LogEdit, ui->progressBar, 100);

    on_PlotButton_clicked();
    Output_File_Write_(ui->LogEdit);
    ui->StartButton->setEnabled(true);
    ui->PlotButton->setEnabled(true);
    ui->Plot2DButton->setEnabled(true);
    ui->ButtonDown->setEnabled(true);
    ui->ButtonUp->setEnabled(true);
    if (!ui->PipedCheckBox->checkState()){
    QFile Param(ProjectPath_+"/ProjectParameters.txt");
    if (Param.open(QFile::WriteOnly|QFile::Truncate) )
    {
               QTextStream cout(&Param);
               if (ui->ExperCheckBox->checkState()){
                   cout << I_p_str << endl
                    << c_e_str << a << " nm" << endl
                    << S_R_str << solvent << " nm" << endl
                    << S_E_D_str << rho_0 << " nm^(-3)" << endl
                    << N_p_str << ui->NphiLabel->text().toInt() << endl
                    << N_z_str << ui->NzLabel->text().toInt() << endl
                    << "----------------------------------------------------------" << endl
                    << C_p_str << endl
                    << H_v_str << Protein.volume.hydrated << " nm^3" << endl
                    << N_C_h_s_str << Protein.square.hydrated << endl
                    << M_v_str << Protein.volume.molecular << " nm^3" << endl
                    << N_c_m_s_str << Protein.square.molecular << endl
                    << C_str << Protein.GetChi() << endl
                    << G_R_str << Protein.RadiusGuinier() <<" nm"<< endl;
               }
               else{
                   cout << I_p_str << endl
//                    << c_e_str << a << " nm" << endl
                    << S_R_str << solvent << " nm" << endl
                    << S_E_D_str << rho_0 << " nm^(-3)" << endl
                    << Q_m_str << Q_min << " nm^(-1)" << endl
                    << Q_M_str << Q_max << " nm^(-1)" << endl
                    << N_Q_str << Q_n << endl
                    << N_p_str << ui->NphiLabel->text().toInt() << endl
                    << N_z_str << ui->NzLabel->text().toInt() << endl
                    << "---------------------------------------------" << endl
                    << C_p_str << endl
                    << H_v_str << Protein.volume.hydrated << " nm^3" << endl
                    << N_C_h_s_str << Protein.square.hydrated << endl
                    << M_v_str << Protein.volume.molecular << " nm^3" << endl
                    << N_c_m_s_str << Protein.square.molecular << endl
                    << G_R_str << Protein.RadiusGuinier() <<" nm"<< endl;
               }
        }
        Param.close();

        if (Param.open(QFile::ReadOnly|QFile::Truncate)){
            ui->ProtInfo->clear();
            QTextStream stream(&Param);
            while (!stream.atEnd()){
                ui->ProtInfo->append(stream.readLine());
            }
        }
        Param.close();
    }
    else{
        QFile Param(ProjectPath_+"/ProjectParameters.txt");
        if (Param.open(QFile::WriteOnly|QFile::Append) )
        {
                   QTextStream cout(&Param);

                       cout << endl
                        << "---------------------------------------------" << endl
                        << "Recalculated using Parallelepipeds file" << endl
                        << I_p_str << endl
                        << Q_m_str << Q_min << " nm^(-1)" << endl
                        << Q_M_str << Q_max << " nm^(-1)" << endl
                        << N_Q_str << Q_n << endl
                        << N_p_str << ui->NphiLabel->text().toInt() << endl
                        << N_z_str << ui->NzLabel->text().toInt() << endl
                        << "---------------------------------------------" << endl
                        << C_p_str << endl
                        << G_R_str << Protein.RadiusGuinier() <<" nm"<< endl;

            }
            Param.close();

            if (Param.open(QFile::ReadOnly|QFile::Truncate)){
                ui->ProtInfo->clear();
                QTextStream stream(&Param);
                while (!stream.atEnd()){
                    ui->ProtInfo->append(stream.readLine());
                }
            }
            Param.close();
        }
    std::time_t Total_TIME = std::time(nullptr)-Start_TIME;
    ui->LogEdit->append("Total Execution Time: "+ QString::number(Total_TIME) +" seconds");
    fclose(stdout);
    ui->statusBar->showMessage("Done", 0);

}

void MainWindow::on_PlotButton_clicked()
{
    QLineSeries *series = new QLineSeries();
    QString way = ProjectPath_+"/curve.txt";
    std::ifstream curve(way.toStdString());
    std::vector<double> temp(2);
    while (!curve.eof()) {
        for (unsigned int i = 0; i<2; ++i)
            curve >> temp[i];
        *series << QPointF(temp[0], temp[1]);
    }
    chart = new QChart();
    chart->legend()->hide();
    chart->addSeries(series);
    chart->setAnimationOptions(QChart::SeriesAnimations);

    QValueAxis *axisX = new QValueAxis();
    axisX->setTitleText("1/nm");
    axisX->setLabelFormat("%.1g");
    axisX->setTickCount(6);
    axisX->setMinorTickCount(4);
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisX);

    QLogValueAxis *axisY = new QLogValueAxis();
    axisY->setTitleText("Intensity");
    axisY->setLabelFormat("%.1g");
    axisY->setBase(10);
    axisY->setMinorTickCount(-1);
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisY);

    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    gridLayout->addWidget(chartView,0,0);
}

void MainWindow::on_actionAbout_triggered()
{
    about = new About (this);
    about->show();
}

void MainWindow::on_ExperCheckBox_clicked()
{
    if (ui->ExperCheckBox->checkState()==false){
        ui->ExperFilePath->setEnabled(false);
        ui->ExperChooseButon->setEnabled(false);
        ui->Qmin->setEnabled(true);
        ui->Qmax->setEnabled(true);
        ui->Qn->setEnabled(true);

    }
    else{
        ui->ExperFilePath->setEnabled(true);
        ui->ExperChooseButon->setEnabled(true);

        ui->Qmin->setEnabled(false);
        ui->Qmax->setEnabled(false);
        ui->Qn->setEnabled(false);
    }
}

void MainWindow::on_ExperChooseButon_clicked()
{
    QString Path = QFileDialog::getOpenFileName(this,
        tr("File of experimental data"), "/", tr("data (*.dat)"));
    ui->ExperFilePath->setText(Path);

    std::ifstream data(Path.toStdString());
    std::vector<double> temp0(3);
    while (!data.eof()) {

        data >> temp0[0] >> temp0[1] >> temp0[2];
        temp0[0] *=10.0;
        Ie_.push_back(temp0);
    }
    ui->LogEdit->append(QString::number(Ie_.size()));
}


void MainWindow::on_Plot2DButton_clicked()
{
    ui->statusBar->showMessage("Calculating Cross Section", 0);
    int x_min = 0, x_max = 0, y_min = 0, y_max = 0;
    QScatterSeries *series = new QScatterSeries();
    series->setMarkerShape(QScatterSeries::MarkerShapeRectangle);
    series->setBorderColor(Qt::blue);
    series->setMarkerSize(2.0);

    QScatterSeries *seriesBorder = new QScatterSeries();
    seriesBorder->setMarkerShape(QScatterSeries::MarkerShapeRectangle);
    seriesBorder->setBorderColor(Qt::red);
    seriesBorder->setMarkerSize(2.0);

    QString wayHyd = ProjectPath_ + "/Hydrated.txt";
    QString wayMol = ProjectPath_ + "/Molecular.txt";

    std::ifstream plot2dHyd(wayHyd.toStdString());
    std::ifstream plot2dMol(wayMol.toStdString());
    
    std::vector<int> temp(4);
    std::vector< std::vector<int> > Protein3DHyd;
    std::vector< std::vector<int> > Protein3DMol;

    int normal_coord = (ui->NormCoordEdit->text()).toInt();

    while (!plot2dHyd.eof()) {
        for (unsigned int i = 0; i < 4; ++i)
            plot2dHyd >> temp[i];
        Protein3DHyd.push_back(temp);
    }

    while (!plot2dMol.eof()) {
        for (unsigned int i = 0; i < 4; ++i)
            plot2dMol >> temp[i];
        Protein3DMol.push_back(temp);
    }


    if (ui->NormalBox->currentText() == "Z"){
        sort(Protein3DHyd.begin(), Protein3DHyd.end(), ZXY_);
        sort(Protein3DMol.begin(), Protein3DMol.end(), ZXY_);

        x_min = Protein3DMol[0][0];
        x_max = Protein3DMol[0][0];
        y_min = Protein3DMol[0][1];
        y_max = Protein3DMol[0][1];

        for (auto &i: Protein3DMol){
            if (normal_coord == i[2])
                *series << QPointF(i[0], i[1]);
            if (i[0] < x_min)
                x_min = i[0];
            if (i[0] > x_max)
                x_max = i[0];

            if (i[1] < y_min)
                y_min = i[1];
            if (i[1] > y_max)
                y_max = i[1];
        }


        x_min = Protein3DHyd[0][0];
        x_max = Protein3DHyd[0][0];
        y_min = Protein3DHyd[0][1];
        y_max = Protein3DHyd[0][1];

        for (auto &i: Protein3DHyd){
            if (normal_coord == i[2] && i[3] == 2)
                *seriesBorder << QPointF(i[0], i[1]);
            if (i[0] < x_min)
                x_min = i[0];
            if (i[0] > x_max)
                x_max = i[0];

            if (i[1] < y_min)
                y_min = i[1];
            if (i[1] > y_max)
                y_max = i[1];
        }

        ui->PositionSlider->setRange(Protein3DHyd[0][2], Protein3DHyd[Protein3DHyd.size()-1][2]);

        chart = new QChart();
        chart->legend()->hide();
        chart->addSeries(series);
        chart->addSeries(seriesBorder);
        chart->setAnimationOptions(QChart::SeriesAnimations);

        QValueAxis *axisX = new QValueAxis();
        axisX->setTitleText("x");
        axisX->setLabelFormat("%i");

        axisX->setRange(x_min-3, x_max+3);

        chart->addAxis(axisX, Qt::AlignBottom);
        series->attachAxis(axisX);
        seriesBorder->attachAxis(axisX);


        QValueAxis *axisY = new QValueAxis();
        axisY->setTitleText("y");
        axisY->setLabelFormat("%i");

        axisY->setRange(y_min-3, y_max+3);

        chart->addAxis(axisY, Qt::AlignLeft);
        series->attachAxis(axisY);
        seriesBorder->attachAxis(axisY);
    }

    else if (ui->NormalBox->currentText() == "Y"){
        sort(Protein3DHyd.begin(), Protein3DHyd.end(), YZX_);
        sort(Protein3DMol.begin(), Protein3DMol.end(), YZX_);

        x_min = Protein3DMol[0][0];
        x_max = Protein3DMol[0][0];
        y_min = Protein3DMol[0][2];
        y_max = Protein3DMol[0][2];

        for (auto &i: Protein3DMol){
            if (normal_coord == i[1])
                *series << QPointF(i[0], i[2]);

            if (i[0] < x_min)
                x_min = i[0];
            if (i[0] > x_max)
                x_max = i[0];

            if (i[2] < y_min)
                y_min = i[2];
            if (i[2] > y_max)
                y_max = i[2];
        }


        x_min = Protein3DHyd[0][0];
        x_max = Protein3DHyd[0][0];
        y_min = Protein3DHyd[0][2];
        y_max = Protein3DHyd[0][2];

        for (auto &i: Protein3DHyd){
            if (normal_coord == i[1] && i[3] == 2)
                *seriesBorder << QPointF(i[0], i[2]);
            if (i[0] < x_min)
                x_min = i[0];
            if (i[0] > x_max)
                x_max = i[0];

            if (i[2] < y_min)
                y_min = i[2];
            if (i[2] > y_max)
                y_max = i[2];
        }
;
        ui->PositionSlider->setRange(Protein3DHyd[0][1], Protein3DHyd[Protein3DHyd.size()-1][1]);


        chart = new QChart();
        chart->legend()->hide();
        chart->addSeries(series);
        chart->addSeries(seriesBorder);
        chart->setAnimationOptions(QChart::SeriesAnimations);

        QValueAxis *axisX = new QValueAxis();
        axisX->setTitleText("x");
        axisX->setLabelFormat("%i");
        axisX->setRange(x_min-3, x_max+3);
        chart->addAxis(axisX, Qt::AlignBottom);
        series->attachAxis(axisX);
        seriesBorder->attachAxis(axisX);

        QValueAxis *axisY = new QValueAxis();
        axisY->setTitleText("z");
        axisY->setLabelFormat("%i");
        axisY->setRange(y_min-3, y_max+3);
        chart->addAxis(axisY, Qt::AlignLeft);
        series->attachAxis(axisY);
        seriesBorder->attachAxis(axisY);
    }

    else if (ui->NormalBox->currentText() == "X"){
        sort(Protein3DHyd.begin(), Protein3DHyd.end(), XYZ_);
        sort(Protein3DMol.begin(), Protein3DMol.end(), XYZ_);

        x_min = Protein3DMol[0][1];
        x_max = Protein3DMol[0][1];
        y_min = Protein3DMol[0][2];
        y_max = Protein3DMol[0][2];

        for (auto &i: Protein3DMol){
            if (normal_coord == i[0])
                *series << QPointF(i[1], i[2]);

            if (i[1] < x_min)
                x_min = i[1];
            if (i[1] > x_max)
                x_max = i[1];

            if (i[2] < y_min)
                y_min = i[2];
            if (i[2] > y_max)
                y_max = i[2];
        }



        x_min = Protein3DHyd[0][1];
        x_max = Protein3DHyd[0][1];
        y_min = Protein3DHyd[0][2];
        y_max = Protein3DHyd[0][2];

        for (auto &i: Protein3DHyd){
            if (normal_coord == i[0] && i[3] == 2)
                *seriesBorder << QPointF(i[1], i[2]);
            if (i[1] < x_min)
                x_min = i[1];
            if (i[1] > x_max)
                x_max = i[1];

            if (i[2] < y_min)
                y_min = i[2];
            if (i[2] > y_max)
                y_max = i[2];
        }

        ui->PositionSlider->setRange(Protein3DHyd[0][0], Protein3DHyd[Protein3DHyd.size()-1][0]);

        chart = new QChart();
        chart->legend()->hide();
        chart->addSeries(series);
        chart->addSeries(seriesBorder);
        chart->setAnimationOptions(QChart::SeriesAnimations);

        QValueAxis *axisX = new QValueAxis();
        axisX->setTitleText("y");
        axisX->setLabelFormat("%i");
        axisX->setRange(x_min-3, x_max+3);
        chart->addAxis(axisX, Qt::AlignBottom);
        series->attachAxis(axisX);
        seriesBorder->attachAxis(axisX);

        QValueAxis *axisY = new QValueAxis();
        axisY->setTitleText("z");
        axisY->setLabelFormat("%i");
        axisY->setRange(y_min-3, y_max+3);
        chart->addAxis(axisY, Qt::AlignLeft);
        series->attachAxis(axisY);
        seriesBorder->attachAxis(axisY);
    }
    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    gridLayout2D->addWidget(chartView,0,0);

    ui->statusBar->showMessage("Cross Section is plotted", 0);
}


bool MainWindow::XYZ_(const std::vector <int>& vec1, const std::vector <int>& vec2) {
    for (int i = 0; i <= int(CubeComponents_::kZ); ++i) {
        if (vec1[i] == vec2[i]) continue;
        return vec1[i] < vec2[i];
    }
    return false;
}

bool MainWindow::YZX_(const std::vector <int>& vec1, const std::vector <int>& vec2) {
    for (int i = 1; i <= int(CubeComponents_::kZ); ++i) {
        if (vec1[i] == vec2[i]) continue;
        return vec1[i] < vec2[i];
    }
    return vec1[int(CubeComponents_::kX)] < vec2[int(CubeComponents_::kX)];
}


bool MainWindow::ZXY_(const std::vector <int>& vec1, const std::vector <int>& vec2) {
    if (vec1[int(CubeComponents_::kZ)] == vec2[int(CubeComponents_::kZ)]) {
        for (int i = 0; i <int(CubeComponents_::kZ); ++i) {
            if (vec1[i] == vec2[i]) continue;
            return vec1[i] < vec2[i];
        }
    }
    return vec1[int(CubeComponents_::kZ)] < vec2[int(CubeComponents_::kZ)];
}

void MainWindow::on_PositionSlider_valueChanged(int value)
{
    ui->NormCoordEdit->setText(QString::number(value));
}

void MainWindow::on_NormCoordEdit_editingFinished()
{
    ui->PositionSlider->setValue( (ui->NormCoordEdit->text()).toInt());
}


void MainWindow::on_ButtonDown_clicked()
{
    int value = (ui->PositionSlider->value());
    ui->PositionSlider->setSliderPosition(--value);
    on_Plot2DButton_clicked();
}

void MainWindow::on_ButtonUp_clicked()
{
    int value = (ui->PositionSlider->value());
    ui->PositionSlider->setSliderPosition(++value);
    on_Plot2DButton_clicked();
}

