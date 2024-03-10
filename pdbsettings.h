#ifndef PDBSETTINGS_H
#define PDBSETTINGS_H

#include <QDialog>
#include <QFile>
#include <QTextStream>
#include <QFileDialog>

namespace Ui {
class PDBSettings;
}

class PDBSettings : public QDialog
{
    Q_OBJECT

public:
    explicit PDBSettings(QWidget *parent = 0);
    ~PDBSettings();

private slots:
    void on_buttonBox_accepted();

    void on_ChooseFolderButton_clicked();

    void on_ChooseFileButton_clicked();

private:
    Ui::PDBSettings *ui;
};

#endif // PDBSETTINGS_H
