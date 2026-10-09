// Persistent workers for independent VBD vertex blocks; no environment-variable configuration.
#ifndef BT_DEFORMABLE_VBD_PARALLEL_H
#define BT_DEFORMABLE_VBD_PARALLEL_H
#include <atomic>
#include <condition_variable>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>
class btDeformableVbdWorkers
{
	std::vector<std::thread> threads;
	std::mutex mutex;
	std::condition_variable ready, complete;
	std::function<void(int)> operation;
	std::atomic<int> next{0};
	int count = 0, generation = 0, remaining = 0;
	bool stop = false;
	void execute()
	{
		for (int begin = next.fetch_add(16); begin < count; begin = next.fetch_add(16))
			for (int i = begin; i < count && i < begin + 16; ++i)
				operation(i);
	}

  public:
	explicit btDeformableVbdWorkers(int workers)
	{
		for (int i = 1; i < workers; ++i)
			threads.emplace_back(
				[this]
				{
					int observed = 0;
					for (;;)
					{
						std::unique_lock<std::mutex> lock(mutex);
						ready.wait(lock, [&] { return stop || generation != observed; });
						if (stop)
							return;
						observed = generation;
						lock.unlock();
						execute();
						lock.lock();
						if (--remaining == 0)
							complete.notify_one();
					}
				});
	}
	~btDeformableVbdWorkers()
	{
		{
			std::lock_guard<std::mutex> lock(mutex);
			stop = true;
		}
		ready.notify_all();
		for (auto &thread : threads)
			thread.join();
	}
	void run(int n, const std::function<void(int)> &fn)
	{
		{
			std::lock_guard<std::mutex> lock(mutex);
			count = n;
			next = 0;
			operation = fn;
			remaining = int(threads.size());
			++generation;
		}
		ready.notify_all();
		execute();
		std::unique_lock<std::mutex> lock(mutex);
		complete.wait(lock, [&] { return remaining == 0; });
		operation = {};
	}
};
inline void btVbdParallelFor(int count, int workers, const std::function<void(int)> &operation)
{
	if (workers <= 1 || count < 64)
	{
		for (int i = 0; i < count; ++i)
			operation(i);
		return;
	}
	static thread_local std::unique_ptr<btDeformableVbdWorkers> pool;
	static thread_local int size = 0;
	if (size != workers)
	{
		pool.reset(new btDeformableVbdWorkers(workers));
		size = workers;
	}
	pool->run(count, operation);
}
// Contiguous batches amortize scheduling for inexpensive geometry operations.
template <class Operation> inline void btVbdParallelGeometry(int count, int workers, const Operation &operation)
{
	const int batch = 256;
	btVbdParallelFor((count + batch - 1) / batch, workers,
					 [&](int block)
					 {
						 const int end = count < (block + 1) * batch ? count : (block + 1) * batch;
						 for (int i = block * batch; i < end; ++i)
							 operation(i);
					 });
}
#endif
