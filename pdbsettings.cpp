#include "pdbsettings.h"
#include "ui_pdbsettings.h"

PDBSettings::PDBSettings(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PDBSettings)
{
    ui->setupUi(this);
    QFile file("Settings.txt");
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    int t = 0;
    QTextStream in(&file);
    QString line= in.readLine();
    while (!line.isNull()) {
        if (t==0) ui->FilePath->setText(line);
        else if (t==1) ui->FolderPath->setText(line);
        line = in.readLine();
        ++t;
    }
    file.close();
}

PDBSettings::~PDBSettings()
{
    delete ui;
}

void PDBSettings::on_buttonBox_accepted()
{
    QString filename="Settings.txt";
    QFile file( filename );
    if (file.open(QFile::WriteOnly|QFile::Truncate) )
    {
        QTextStream cout( &file );
       cout<<ui->FilePath->text()<<endl;
       cout<<ui->FolderPath->text()<<endl;
    }
    file.close();
    QDialog::reject();
}

void PDBSettings::on_ChooseFolderButton_clicked()
{
   QString Path = QFileDialog::getExistingDirectory(this);
   ui->FolderPath->setText(Path);
}

void PDBSettings::on_ChooseFileButton_clicked()
{
    ui->FilePath->setText(QFileDialog::getOpenFileName(this,
              tr("Open the file, which contains descriptions of atoms groups:"), "", tr("Text (*.txt)")));
}
