package mairaDatabase.utils;

import java.util.Collections;
import java.util.List;

public class Statistics {

	public static double getMean(List<Double> data) {
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

	public static double getVariance(List<Double> data) {
		double mean = getMean(data);
		double temp = 0;
		for (double a : data)
			temp += (a - mean) * (a - mean);
		return temp / (data.size() - 1);
	}

	public static double getStdDev(List<Double> data) {
		return Math.sqrt(getVariance(data));
	}

	public static double getMedian(List<Double> data) {
		Collections.sort(data);
		if (data.size() % 2 == 0)
			return (data.get((data.size() / 2) - 1) + data.get(data.size() / 2)) / 2.0;
		return data.get(data.size() / 2);
	}

	public static double getMin(List<Double> data) {
		Collections.sort(data);
		return data.get(0);
	}

	public static double getMax(List<Double> data) {
		Collections.sort(data);
		return data.get(data.size() - 1);
	}

}
