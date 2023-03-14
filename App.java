import java.awt.geom.Point2D;
import java.io.*;
import java.util.*;

import jimena.binaryrn.RegulatoryNetwork;
import jimena.perturbation.OnOffPerturbation;
import jimena.simulationmethods.NormalizedHillCubeInterpolationMethod;
import jimena.simulationmethods.SQUADInterpolationMethod;
/**
 * 
 */
public class App {
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        String PARAMETER_INPUT = "jimena-app/ParameterInputs.csv";
        
        RegulatoryNetwork network = new RegulatoryNetwork();
        try {
            String fileName = args[0];
            network.loadYEdFile(new File(fileName));
        } catch (Exception e) {
            e.printStackTrace();
        }

        String [] nodeNames = network.getNodeNames();
        deleteWhiteSpacesFromNodeNames(nodeNames);
        ArrayList<String> headers = new ArrayList<String>(Arrays.asList(nodeNames));
        headers.add(0,"time");
        

        BufferedReader csvReader = new BufferedReader(new FileReader(PARAMETER_INPUT));
        String row;
        Integer method = 0;
        Double dt = null;
        Double maxt = null;
        Double maxSpeed = null;
        Double minSimulationTimeBtwLogs = null;
        Double minTimeBtwNotification = null;
        int currentLine = 0;
        while ((row = csvReader.readLine()) != null) {
            if(currentLine == 0) {
                // Read params for simulation and store them.
                String[] simulationParameters = row.split(",");
                method = Integer.parseInt(simulationParameters[0]);
                dt = Double.parseDouble(simulationParameters[1]);
                maxt = Double.parseDouble(simulationParameters[2]);
                maxSpeed = Double.parseDouble(simulationParameters[3]);
                minSimulationTimeBtwLogs = Double.parseDouble(simulationParameters[4]);
                minTimeBtwNotification  = Double.parseDouble(simulationParameters[5]);

                currentLine++;
            } else {
                // Read perturbationParams from file and add network
                String[] perturbationParams = row.split(",");
                int index = Integer.parseInt(perturbationParams[0]);
                int nodeIndex = Integer.parseInt(perturbationParams[1]);
                int start = Integer.parseInt(perturbationParams[2]);
                int end = Integer.parseInt(perturbationParams[3]);
                double value = Double.parseDouble(perturbationParams[4]);
                addPerturbation(network, index, nodeIndex, start, end, value);
                currentLine++;
            }
        }
        csvReader.close();

        simulateNetworkWithParams(network, method, dt, maxt, maxSpeed, minSimulationTimeBtwLogs, minTimeBtwNotification);

        // Create 2d array inorder to export logs as CSV
        int columnSize = headers.size();
        int rowSize = network.getNetworkNodes()[0].getLog().size();
        String[][] csvTable = new String[rowSize][columnSize];

        // Write the logs to the 2d array
        for(int i = 0; i < nodeNames.length; i++) {
            int counter = 0;
            for (Point2D.Double logEntry : network.getNetworkNodes()[i].getLog()) {
                csvTable[counter][0] = Double.toString(logEntry.getX());
                csvTable[counter][i + 1] = Double.toString(logEntry.getY());
                counter++;
            }
        }

        BufferedWriter writer = new BufferedWriter(new FileWriter("jimena_time_series_data.csv", true));

        // Write headers to csv file
        writer.append(concatStringArray(headers,","));
        writer.append("\n");
        // Write logs to csv file
        for(int i = 0; i < csvTable.length; i++) {
            System.out.println("row: " + concatStringArray(csvTable[i], ","));
            writer.append(concatStringArray(csvTable[i], ","));
            writer.append("\n");
        }
        
        writer.close();
            System.exit(0);
    }



    /** Simulate network with params
     *
     * @param network Network object
     * @param method Perturbation number
     * @param dt Index of node
     * @param maxt Perturbation param
     * @param maxSpeed Perturbation param
     * @param minSimulationTimeBtwLogs Perturbation param
     * @param minTimeBtwNotifications
     *
     */
    private static void simulateNetworkWithParams(RegulatoryNetwork network, int method, double dt, double maxt, double maxSpeed, double minSimulationTimeBtwLogs, double minTimeBtwNotifications) {
        network.simulate(new SQUADInterpolationMethod(),
                dt,
                maxt,
                maxSpeed,
                minSimulationTimeBtwLogs,
                minTimeBtwNotifications,
                null
        );
    }
    /** Add selected perturbation to the network node
     *
     * @param network Network object
     * @param index Perturbation number
     * @param nodeIndex Index of node
     * @param start Perturbation param
     * @param end Perturbation param
     * @param value Perturbation param
     *
     */
    private static void addPerturbation(RegulatoryNetwork network, int index, int nodeIndex, int start, int end, double value) {
        if(index == 0) {
            network.getNetworkNodes()[nodeIndex].getPerturbations().add(new OnOffPerturbation(start, end, value));
            
        }
            
    }

    /**
     * Convert string array to strings
     */
    public static String concatStringArray(String[] arr, String separator) {
        if (null == arr || 0 == arr.length) return "";

        String row = "";

        for (int i = 0; i < arr.length; i++) {
            if(i == arr.length - 1) {
                row += arr[i];
            } else {
                row += arr[i] + separator;
            }
        }
        return row;
    }

    /**
     * Convert string arraylist to strings
     */
    public static String concatStringArray(ArrayList<String> arr, String separator) {
        if (null == arr || 0 == arr.size()) return "";

        String row = "";

        for (int i = 0; i < arr.size(); i++) {
            if(i == arr.size() - 1) {
                row += arr.get(i);
            } else {
                row += arr.get(i) + separator;
            }
        }
        System.out.println(row);
        return row;
    }

    public static void deleteWhiteSpacesFromNodeNames(String[] arr) {
        for(int i = 0; i < arr.length; i++) {
            arr[i] = arr[i].replaceAll("\\s", ""); // using built in method
        }
    }
}

