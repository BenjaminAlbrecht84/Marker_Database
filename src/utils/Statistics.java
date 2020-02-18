package utils;

import java.util.ArrayList;
import java.util.Collections;

public class Statistics {

	public static double getMean(ArrayList<Double> data) {
		double sum = 0.0;
		for (double a : data)
			sum += a;
		return sum / data.size();
	}

	public static double getMean(long[] vector) {
		double sum = 0.0;
		for (long a : vector)
			sum += a;
		return sum / vector.length;
	}

	public static double getVariance(ArrayList<Double> data) {
		double mean = getMean(data);
		double temp = 0;
		for (double a : data)
			temp += (a - mean) * (a - mean);
		return temp / (data.size() - 1);
	}

	public static double getStdDev(ArrayList<Double> data) {
		return Math.sqrt(getVariance(data));
	}

	public static double getMedian(ArrayList<Double> data) {
		Collections.sort(data);
		if (data.size() % 2 == 0)
			return (data.get((data.size() / 2) - 1) + data.get(data.size() / 2)) / 2.0;
		return data.get(data.size() / 2);
	}

	public static double getMin(ArrayList<Double> data) {
		Collections.sort(data);
		return data.get(0);
	}

	public static double getMax(ArrayList<Double> data) {
		Collections.sort(data);
		return data.get(data.size() - 1);
	}

}
