package mairaDatabase.utils;

import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicLong;

public class ResourceLoader {
	
	private long maxProgress, lastProgress = 0;
	private AtomicLong progress = new AtomicLong();
	private long time = System.currentTimeMillis();

	private CountDownLatch latch;
	private ExecutorService executor;
	
	public void runThreads(int size, List<Runnable> threads, long maxProgress) {
		setMaxProgress(maxProgress);
        setSize(size);
        setTime();
        runInParallel(threads);
        shutdown();
        reportFinish();
        System.out.println("Runtime: " + getUptime());
	}
	
	public void setSize(int size) {
		executor = Executors.newFixedThreadPool(size);
	}
	
	public void setTime() {
		time = System.currentTimeMillis();
	}
	
	public void setMaxProgress(long maxProgress) {
		this.maxProgress = maxProgress;
	}
	
	public void shutdown() {
		if(executor != null)
			executor.shutdown();
	}
	
	public void countDown() {
		if(latch != null)
			latch.countDown();
	}
	
	public void runInParallel(List<Runnable> threads) {
		latch = new CountDownLatch(threads.size());
		for (Runnable t : threads)
			executor.submit(t);
		try {
			latch.await();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}

	public synchronized void reportProgress(long delta) {
		progress.getAndAdd(delta);
		int p = ((int) ((((double) progress.get() / (double) maxProgress)) * 100) / 5) * 5;
		if (p > lastProgress && p < 100) {
			lastProgress = p;
			System.out.print(p + "% (" + getUptime() + ") ");
		}
	}

	public void reportFinish() {
		progress.set(0);
		lastProgress = 0;
		System.out.print(100 + "%\n");
	}

	public String getUptime() {
		long runtime = (System.currentTimeMillis() - time) / 1000;
		return runtime + "s";
	}
	
	public void reportRuntime() {
		System.out.println("Runtime: " + getUptime());
	}

	public CountDownLatch getLatch() {
		return latch;
	}
	
}
