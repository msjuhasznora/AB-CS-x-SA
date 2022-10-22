// az antibiotikum túl gyors, túl erős -- realisztikusabbá tenni, sztochasztikusan marad pár egy darabig?
// efficacy-k
// késleltetett immunválasz megoldása (2-3 nap?) --> ne instant reagáljon
// immunválasz lassú lecsengése

// mit szeretnénk kihozni:
// 1) ha 2-3 napon belül kap antibiotikumot, és mást nem, akkor nem károsodik a porc jelentősen
// 2) ha az 5-6. napon jut csak antibiotikumhoz, és máshoz nem, akkor jelentősen károsodik a porc
// 3) ha az 5.-6. napon jut csak antibiotikumhoz, viszont ekkor kortikoszteroidot is kap, akkor nem károsodik a porc

package CorticosteroidsAntibiotics;

import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import HAL.Tools.FileIO;
import HAL.Rand;
import static HAL.Util.*;

import java.io.File;

public class CorticosteroidsAntibiotics{

    public static final int BIG_VALUE = (Integer.MAX_VALUE-1) / 2;

    public static void main(String[] args) {

        String singularOrSweep = "singular"; // "singular" or "sweep"
        String inVivoOrInVitro = "inVivo";

        int y = 100, x = 100, visScale = 3;
        boolean isAntibiotics = true;
        boolean isCorticosteroids = false;

        GridWindow win = new GridWindow("Cellular state space, bacterial concentration, immune reaction.", x*3, y, visScale,true);

        if (singularOrSweep.equals("singular")){

            NewExperiment experiment = new NewExperiment(x, y, visScale, new Rand(1), isAntibiotics, isCorticosteroids, 1*2*12*60, 0.2, inVivoOrInVitro, 110.0);
            experiment.numberOfTicks = experiment.numberOfTicksDelay + experiment.numberOfTicksDrug;

            experiment.Init();
            double remainingHealthyCells = experiment.RunExperiment(win);

            if (isAntibiotics) {
                System.out.println(inVivoOrInVitro.equals("inVivo") ? "In vivo drug source [ng / ml]: " + experiment.antibiotics.drugSourceStomach : "In vitro drug concentration [nM]: " + experiment.antibiotics.NgPerMlToNanomolars(experiment.antibiotics.inVitroAntibioticsCon));
            }
            System.out.println("Remaining healthy cells: " + remainingHealthyCells);

        } else if (singularOrSweep.equals("sweep")) {

            String collectiveOutputDir = collectiveOutputDir(isAntibiotics, isCorticosteroids);
            FileIO collectiveOutFile = new FileIO(collectiveOutputDir.concat("/").concat("collectiveRemainingHealthyCells").concat(".csv"), "w");
            String collectiveResults;

            for (double staphyloDiffCoeffSweep = 0.00625; staphyloDiffCoeffSweep < 50; staphyloDiffCoeffSweep *= 2) {

                collectiveResults = "";
                collectiveResults += staphyloDiffCoeffSweep + ", ";

                for (double damageRateSweep = 0.0; damageRateSweep < 100.0; damageRateSweep += 5.0) {

                    NewExperiment experiment = new NewExperiment(x, y, visScale, new Rand(1), isAntibiotics, isCorticosteroids, BIG_VALUE, staphyloDiffCoeffSweep, inVivoOrInVitro, damageRateSweep);
                    experiment.numberOfTicks = experiment.numberOfTicksDelay + experiment.numberOfTicksDrug;

                    experiment.Init();
                    double remainingHealthyCells = experiment.RunExperiment(win);

                    collectiveResults += remainingHealthyCells + ", ";

                }
                collectiveOutFile.Write(collectiveResults + "\n");
                System.out.println(collectiveResults);
            }

            collectiveOutFile.Close();

        } else {
            System.out.println("The only options are singular and sweep.");
        }

        win.Close();

    }

    public static String collectiveOutputDir(boolean isAntibiotics, boolean isCorticosteroid){

        java.util.Date now = new java.util.Date();
        java.text.SimpleDateFormat dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
        String date_time = dateFormat.format(now);
        String projPath = PWD() + "/output/StaphyloExperiments";
        if (isAntibiotics == false){
            projPath += "/noDrug";
        } else if (isCorticosteroid == true){
            projPath += "/steroidBoostedAB";
        } else {
            projPath += "/AntibioticsOnly";
        }
        String collectiveOutputDir = projPath + "/" + date_time + "__collective";
        new File(collectiveOutputDir).mkdirs();

        return collectiveOutputDir;

    }

}

class CorticosteroidDrug{

    // general properties
    double EC50 = 62; // in nM = nanoMolars, [nM] = 10^-9 [mol/L]; https://www.fda.gov/media/155050/download
    double molarMassDrug = 499.535;

    // in vitro properties
    double inVitroAntibioticsCon = 5;

    // in vitro constructor
    public CorticosteroidDrug(double inVitroAntibioticsCon, boolean isCorticosteroids){

        if (isCorticosteroids == true){
            this.inVitroAntibioticsCon = inVitroAntibioticsCon;
        } else {
            this.inVitroAntibioticsCon = 0.0;
        }

    }

    // in vivo properties

    double drugDecay = 0.013; // see the paxlovid_config.nb Mathematica notebook
    double drugSourceStomach = 1800; // see the paxlovid_config.nb Mathematica notebook
    double drugDecayStomach = 0.015; // see the paxlovid_config.nb Mathematica notebook

    // in vivo constructor
    public CorticosteroidDrug(){

    }

    double CorticosteroidEfficacy(double steroidNow){

        // double drugVirusProdEff = 7000 * Math.pow(drugNow, 2)/(1+7000*Math.pow(drugNow,2));
        // drugNow is the drug concentration in nanograms / ml
        // drugNow needs to be converted to [nM]s, as IC50 is given in [nM]s
        double steroidNowInNanoMolars = NgPerMlToNanomolars(steroidNow);
        double steroidEfficacy = 1 / ( 1 + (EC50 / StochasticDrug(steroidNowInNanoMolars)));

        return steroidEfficacy;

    }

    double StochasticDrug(double drug){

        double stdDevOfGaussian = drug / 100;

        Rand random = new Rand();
        double stochasticDrug = random.Gaussian(drug, stdDevOfGaussian);
        stochasticDrug = stochasticDrug > 0.0 ? stochasticDrug : 0.0;

        return stochasticDrug;

    }

    double NgPerMlToNanomolars(double drugNow){
        return drugNow * Math.pow(10,3) / molarMassDrug;
    }

}

class AntibioticsDrug {

    // general properties
    double EC50 = 62; // in nM = nanoMolars, [nM] = 10^-9 [mol/L]; https://www.fda.gov/media/155050/download
    double molarMassDrug = 499.535;

    // in vitro properties
    double inVitroAntibioticsCon = 5;

    // in vitro constructor
    public AntibioticsDrug(double inVitroAntibioticsCon, boolean isAntibiotics){

        if (isAntibiotics == true){
            this.inVitroAntibioticsCon = inVitroAntibioticsCon;
        } else {
            this.inVitroAntibioticsCon = 0.0;
        }

    }

    // in vivo properties

    double drugDecay = 0.013; // see the paxlovid_config.nb Mathematica notebook
    double drugSourceStomach = 1800; // see the paxlovid_config.nb Mathematica notebook
    double drugDecayStomach = 0.015; // see the paxlovid_config.nb Mathematica notebook

    // in vivo constructor
    public AntibioticsDrug(){

    }

    double AntibioticsEfficacy(double antibioticsNow){

        // double drugVirusProdEff = 7000 * Math.pow(drugNow, 2)/(1+7000*Math.pow(drugNow,2));
        // drugNow is the drug concentration in nanograms / ml
        // drugNow needs to be converted to [nM]s, as IC50 is given in [nM]s
        double drugNowInNanoMolars = NgPerMlToNanomolars(antibioticsNow);
        double antibioticsEfficacy = 1 / ( 1 + (EC50 / StochasticDrug(drugNowInNanoMolars)));

        return antibioticsEfficacy;

    }

    double StochasticDrug(double drug){

        double stdDevOfGaussian = drug / 100;

        Rand random = new Rand();
        double stochasticDrug = random.Gaussian(drug, stdDevOfGaussian);
        stochasticDrug = stochasticDrug > 0.0 ? stochasticDrug : 0.0;

        return stochasticDrug;

    }

    double NgPerMlToNanomolars(double drugNow){
        return drugNow * Math.pow(10,3) / molarMassDrug;
    }

}

class NewExperiment extends AgentGrid2D<Cells>{

    public int x = 100;
    public int y = 100;
    public int visScale = 3;
    public int numberOfTicksDelay;
    public int numberOfTicksDrug;
    public int numberOfTicks;
    public PDEGrid2D bacterialCon;
    public PDEGrid2D immuneResponseLevel; // similar to T-cell concentrations, but more generic
    public double antibioticsCon = 0;
    public double antibioticsConStomach = 0;
    public double corticosteroidCon = 0;
    public double corticosteroidConStomach = 0;
    public Rand rn;
    public double[] cellularBacterialCon = new double[length];
    public double[] cellularImmuneResponseLevel = new double[length];

    public double fixedDamageRate;

    public double stapyloReproductionRate = Math.pow(10,-3);
    public double damageRate = 1.01 * Math.pow(10,-7); // beta in the ODE
    public double deathProb = 7.02 * Math.pow(10,-4); // P_D
    public double staphyloDiffCoeff = 0.2; // D_V [sigma^2 / min]

    AntibioticsDrug antibiotics;
    CorticosteroidDrug corticosteroid;

    String inVivoOrInVitro = "inVivo";

    public double immuneResponseDecay = 0.0005;
    public double immuneResponseDiffCoeff = 0.1;

    public boolean isAntibiotics = true;
    public boolean isCorticosteroid = true;

    public FileIO outFile;
    public FileIO paramFile;
    public FileIO concentrationsFile;
    public String outputDir;

    public NewExperiment(int xDim, int yDim, int visScale, Rand rn, boolean isAntibiotics, boolean isCorticosteroid, int numberOfTicksDelay, double staphyloDiffCoeff, String inVivoOrInVitro, double fixedDamageRate){

        super(xDim, yDim, Cells.class);
        this.x = xDim;
        this.y = yDim;
        this.visScale = visScale;
        this.numberOfTicksDelay = numberOfTicksDelay;
        this.staphyloDiffCoeff = staphyloDiffCoeff;
        this.rn = rn;

        this.inVivoOrInVitro = inVivoOrInVitro;
        this.isAntibiotics = isAntibiotics;
        this.isCorticosteroid = isCorticosteroid;

        this.fixedDamageRate = fixedDamageRate;

        if (inVivoOrInVitro.equals("inVivo")) {

            this.antibiotics = new AntibioticsDrug();
            this.corticosteroid = new CorticosteroidDrug();
            this.numberOfTicksDrug = 10 * 24 * 60; // we administer antibiotics for 5 days, i.e. 5*24*60 minutes

        } else {

            this.antibiotics = new AntibioticsDrug(5.0, isAntibiotics);
            this.corticosteroid = new CorticosteroidDrug(5.0, isCorticosteroid);
            this.numberOfTicksDrug = 4 * 24 * 60; // 4 days

        }

        bacterialCon = new PDEGrid2D(xDim, yDim);
        immuneResponseLevel = new PDEGrid2D(xDim, yDim);
        bacterialCon.Update();
        immuneResponseLevel.Update();

        outputDir = this.OutputDirectory();
        outFile = new FileIO(outputDir.concat("/").concat("Out").concat(".csv"), "w");
        paramFile = new FileIO(outputDir.concat("/").concat("Param").concat(".csv"), "w");
        concentrationsFile = new FileIO(outputDir.concat("/").concat("concentrations").concat(".csv"), "w");

    }

    public void Init(){

        WriteHeader();

        int initialPlace = length / 2 - length % 2 + x / 2 - x % 2;

        for (int i = 0; i < length; i++){
            Cells c = NewAgentSQ(i);
            c.CellInit(true,false, false);
            if (i == initialPlace){
                bacterialCon.Add(c.Isq(), 10.0);
                bacterialCon.Update();
            }
        }
    }

    public double RunExperiment(GridWindow win){

        double[] cellCounts = CountCells();
        // System.out.println(cellCounts[0]+", " + cellCounts[1] + ", " + cellCounts[2]);

        for (int tick = 0; tick < this.numberOfTicks; tick ++){

            if ((numberOfTicksDelay == CorticosteroidsAntibiotics.BIG_VALUE) && (fixedDamageRate <= 100) && (((cellCounts[1] + cellCounts[2]) / 40000) * 100 >= fixedDamageRate)){
                this.numberOfTicksDelay = tick - 1;
                this.numberOfTicks = this.numberOfTicksDelay + this.numberOfTicksDrug;
                System.out.println("Diffusion coeff.: " + this.staphyloDiffCoeff + ". Damage info: " + fixedDamageRate + " percent damage was found at tick " + numberOfTicksDelay + ".");
            }

            TimeStep(tick);
            DrawModel(win);

            if( tick > 0 && ( (tick % (24*60)) == 0 ))
                win.ToPNG(outputDir + "day" + Integer.toString(tick/(24*60)) + ".jpg");

            double totalBacterialCon = TotalBacterialCon();
            double totalImmuneResponseLevel = TotalImmuneResponseLevel();
            cellCounts = CountCells();
            concentrationsFile.Write(totalBacterialCon + "," + totalImmuneResponseLevel + "," + antibioticsCon + "," + corticosteroidCon + "\n");
            outFile.Write(tick +"," + cellCounts[0] + "," + cellCounts[1]+
                    "," + cellCounts[2] + "," + totalBacterialCon + "," + totalImmuneResponseLevel + "," + antibioticsCon + "," + corticosteroidCon + "\n");

        }

        CloseFiles();

        return cellCounts[0];

    }

    double[] CountCells(){

        double healthyCells = 0, damagedCells = 0, deadCells = 0;
        double[] cellCount = new double[3];
        for (Cells cell: this){
            if (cell.CellType == 0){
                healthyCells += 1;
            }
            else if (cell.CellType == 1 ){
                damagedCells += 1;
            }
            else if (cell.CellType == 2){
                deadCells += 1;
            }
        }

        cellCount[0] = healthyCells;
        cellCount[1] = damagedCells;
        cellCount[2] = deadCells;

        return cellCount;
    }

    void TimeStepImmune(int tick){

        // decay of the immuneResponseLevel
        for (Cells cell : this){
            double immuneResponseNow = immuneResponseLevel.Get(cell.Isq());
            immuneResponseLevel.Add(cell.Isq(), -immuneResponseDecay * immuneResponseNow);
        }
        immuneResponseLevel.Update();

        // immune response level increases where bacterial concentration is high
        for (Cells cell : this){
            double addedImmuneResponseLevel = ImmuneResponseSource(tick, cell);
            double currentImmuneResponseLevel = immuneResponseLevel.Get(cell.Isq());
            double newImmuneResponseLevel = addedImmuneResponseLevel + currentImmuneResponseLevel;
            immuneResponseLevel.Set(cell.Isq(), newImmuneResponseLevel);
        }

        immuneResponseLevel.DiffusionADI(immuneResponseDiffCoeff);
        immuneResponseLevel.Update();

    }

    void TimeStepStaphylo(int tick){

        // bacterial reproduction
        for (Cells cell : this){
            bacterialCon.Add(cell.Isq(), stapyloReproductionRate * bacterialCon.Get(cell.Isq()));
        }

        bacterialCon.DiffusionADI(staphyloDiffCoeff);
        bacterialCon.Update();

        // removal of bacteria
        for (Cells cell : this){
            // double removalEfficacy = 2/(1+Math.exp(100*drugNow));
            // double removalEfficacy = 100*Math.pow(drugNow, 2)/(1+100*Math.pow(drugNow,2));
            double drugBacterialRemovalEff = 1 / (1 + 1/(Math.pow(antibioticsCon / 100, 2)));
            double immuneBacterialRemovalEff = 2 * 1 / (1 + 1/(Math.pow(immuneResponseLevel.Get(cell.Isq()),2)));
            bacterialCon.Add(cell.Isq(), -drugBacterialRemovalEff * bacterialCon.Get(cell.Isq()));
            bacterialCon.Add(cell.Isq(), -immuneBacterialRemovalEff * bacterialCon.Get(cell.Isq()));
        }
        bacterialCon.Update();

    }

    void TimeStepDrug(int tick){

        if (this.inVivoOrInVitro.equals("inVitro")){

            this.antibioticsCon = this.antibiotics.inVitroAntibioticsCon;
            this.corticosteroidCon = this.corticosteroid.inVitroAntibioticsCon;

        } else if (this.inVivoOrInVitro.equals("inVivo")){

            // decay of the drug
            this.antibioticsCon -= this.antibiotics.drugDecay * this.antibioticsCon;
            this.corticosteroidCon -= this.corticosteroid.drugDecay * this.corticosteroidCon;

            // decay of the drug in the stomach
            // and appearance of the drug at the lung epithelial cells
            double transferQuantity = this.antibiotics.drugDecayStomach * this.antibioticsConStomach;
            this.antibioticsCon += transferQuantity;
            this.antibioticsConStomach -= transferQuantity;

            transferQuantity = this.corticosteroid.drugDecayStomach * this.corticosteroidConStomach;
            this.corticosteroidCon += transferQuantity;
            this.corticosteroidConStomach -= transferQuantity;

            // drug appearance in the stomach
            this.antibioticsConStomach += AntibioticsSourceStomach(tick);
            this.corticosteroidConStomach += CorticosteroidSourceStomach(tick);

        } else {
            System.out.println("inVitro and inVivo are the only two choices currently.");
        }

    }

    void TimeStepCells(int tick){

        for (Cells cell : this){
            cell.CellDamage();
        }
        for (Cells cell : this) {
            cell.CellDeath();
        }

    }

    void TimeStep(int tick){

        TimeStepImmune(tick);
        TimeStepDrug(tick); // both antibiotics and steroid
        TimeStepStaphylo(tick);
        TimeStepCells(tick);

    }

    double TotalBacterialCon() {

        double totalBacterialCon = 0;
        for (int i = 0; i < length; i++){
            cellularBacterialCon[i] = bacterialCon.Get(i);
        }
        for (double bacterialConInCell : cellularBacterialCon ){
            totalBacterialCon = totalBacterialCon + bacterialConInCell;
        }
        return totalBacterialCon;
    }

    double TotalImmuneResponseLevel() {

        double totalImmuneResponseLevel = 0;
        for (int i = 0; i < length; i++){
            cellularImmuneResponseLevel[i] = immuneResponseLevel.Get(i);
        }
        for (double immuneResponseInCell : cellularImmuneResponseLevel ){
            totalImmuneResponseLevel = totalImmuneResponseLevel + immuneResponseInCell;
        }
        return totalImmuneResponseLevel;
    }

    double ImmuneResponseSource(int tick, Cells cell){

        // todo: improve
        double currentBacterialCon = bacterialCon.Get(cell.Isq());
        return 1 / (1 + 1/(Math.pow(currentBacterialCon,2)));

    }

    double AntibioticsSourceStomach(int tick){

        if ((tick > numberOfTicksDelay) && (isAntibiotics == true) && (((tick - numberOfTicksDelay) % (12 * 60)) == 1)) {
            return this.antibiotics.drugSourceStomach;
        } else {
            return 0.0;
        }
    }

    double CorticosteroidSourceStomach(int tick){

        if ((tick > numberOfTicksDelay) && (isCorticosteroid == true) && (((tick - numberOfTicksDelay) % (12 * 60)) == 1)) {
            return this.corticosteroid.drugSourceStomach;
        } else {
            return 0.0;
        }
    }

    void WriteHeader(){

        outFile.Write("tick, healthy cells, damaged cells, dead cells, "
                + "total bacterial conc., total immune response, total AB, total CS \n");
    }

    void CloseFiles(){

        outFile.Close();
        paramFile.Close();
        concentrationsFile.Close();

    }

    String OutputDirectory(){

        java.util.Date now = new java.util.Date();
        java.text.SimpleDateFormat dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
        String date_time = dateFormat.format(now);
        String projPath = PWD() + "/output/StaphyloExperiments";
        if (this.isAntibiotics == false){
            projPath += "/noDrug";
        } else if (this.isCorticosteroid == true){
            projPath += "/steroidBoostedAB";
        } else {
            projPath += "/AntibioticsOnly";
        }

        double drugInfo;
        if (this.isAntibiotics == false){
            drugInfo = 0.0;
        } else if (this.inVivoOrInVitro.equals("inVitro")) {
            drugInfo = this.antibiotics.NgPerMlToNanomolars(this.antibiotics.inVitroAntibioticsCon);
        } else {
            drugInfo = this.antibiotics.drugSourceStomach;
        }

        String outputDir = projPath + "/" + date_time + this.inVivoOrInVitro + drugInfo + "__diff" + this.staphyloDiffCoeff;
        if(this.fixedDamageRate < 100.0){
            outputDir +=  "__damagerate" + this.fixedDamageRate + "/";
        } else {
            outputDir +=  "__delaytime" + this.numberOfTicksDelay + "/";
        }

        new File(outputDir).mkdirs();

        return outputDir;

    }

    void DrawModel(GridWindow vis){

        for (int i = 0; i < length; i++) {
            Cells drawMe = GetAgent(i);

            if (drawMe == null){
                vis.SetPix(i, RGB256(255, 255, 255));
            }
            else{
                if (drawMe.CellType == 0){		// Healthy cells
                    vis.SetPix(i, RGB256(119, 198, 110));
                }
                else if (drawMe.CellType == 1){   // Infected cells
                    vis.SetPix(i, RGB256(124, 65, 120));
                }
                else if (drawMe.CellType == 2){
                    vis.SetPix(i, RGB(0, 0, 0));
                }
                else if (drawMe.CellType == 3){
                    vis.SetPix(i, RGB(255, 255, 255));
                }
            }

            vis.SetPix(ItoX(i) + xDim, ItoY(i), HeatMapGRB(10*bacterialCon.Get(i)));

            vis.SetPix(ItoX(i) + 2 * xDim, ItoY(i), HeatMapRBG(10*immuneResponseLevel.Get(i)));
        }
    }

}

class Cells extends AgentSQ2Dunstackable<NewExperiment>{
    int CellType;

    public void CellInit(boolean isHealthy, boolean isDamaged,
                         boolean isDead){

        if(isHealthy == true){
            this.CellType = 0;
        }
        else if(isDamaged == true){
            this.CellType = 1;
        }
        else if(isDead == true) {
            this.CellType = 2;
        }
    }

    public void CellDamage(){

        // nor antibiotics nor corticosteroids help the cells directly
        // difference w.r.t. certain antivirals, which prevent the cells from getting infected

        // chondral damage: depends on the immune response strength (a strong immune response kills bacteria, but damages the cartilage as well)

        double immuneConAtCell = G.immuneResponseLevel.Get(Isq());

        double damageProb = G.damageRate * G.xDim * G.yDim;
        double effectiveDamageProb = damageProb * immuneConAtCell;

        if (this.CellType == 0){ // healthy cell
            if (G.rn.Double() < effectiveDamageProb) {
                this.CellType = 1;
            }
        }
    }

    public void CellDeath(){

        if (this.CellType == 1) { // damaged
            if(G.rn.Double() < G.deathProb){
                this.CellType = 2;
            }
        }
    }
}
