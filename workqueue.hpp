#include <functional>
#include <future>
#include <deque>
#include <vector>
#include <cassert>
#include <thread>
#include <iostream>

class WorkQueue {
public:
    // Initialize a workqueue with a worklist and a requested number of workers.
    //
    // We own the memory for the worklist after this.
    //
    // If numWorkers is less than 0, use the number of CPU threads available on the
    // platform.
    //
    explicit WorkQueue(std::vector<std::function<void()>> &worklist, int numWorkers = -1) :
        m_work(worklist)
    {
        if (numWorkers < 1) {
            numWorkers = std::thread::hardware_concurrency() + 1;
        }

        int itemsPerWorker = worklist.size() / numWorkers;
        int leftover = worklist.size() % numWorkers;

        int offset = 0;
        int start = 0;
        int end = 0;
        
        for (int i = 0; i < numWorkers; ++i) {
            offset = (int)(leftover > i);
            end = start + (itemsPerWorker - 1) + offset;
            m_workers.emplace_back(std::thread(&WorkQueue::doWork, this, i, start, end));
            start = end + 1;
        }

    }

    // Stop processing work right away and dispose of threads
    virtual ~WorkQueue() {
        abort();
    }

    // Stop processing work right away and dispose of threads
    void abort() {
        m_exit = true;
        m_finish_work = false;
        m_signal.notify_all();

        joinAll();

        {
            std::lock_guard<std::mutex> lg(m_mutex);
            m_work.clear();
        }
    }

    // Finish all work and then dispose of threads afterwards
    void waitForCompletion() {
        m_exit = true;
        m_finish_work = true;
        m_signal.notify_all();

        joinAll();
    }

private:
    std::vector<std::function<void()>> m_work;
    std::vector<std::thread> m_workers;
    std::mutex m_mutex;
    std::condition_variable m_signal;
    std::atomic<bool> m_exit{ false };
    std::atomic<bool> m_finish_work{ true };  // override m_exit until the work is done

    // Thread main loop
    void doWork(int thr_id, int startIdx, int endIdx) {
        std::unique_lock<std::mutex> ul(m_mutex); // constructed locked
        std::cout << "Th[" << thr_id << "]: " << startIdx << "-" << endIdx << std::endl;

        int curItem = startIdx;
        while (!m_exit || (m_finish_work && (curItem <= endIdx))) {
            if (curItem <= endIdx) {
                std::function<void()> work(std::move(m_work[curItem]));
                ul.unlock();
                work();
                ul.lock();
                ++curItem;
            }
            else{
                m_signal.wait(ul);
            }
        }
    }

    void joinAll() {
        for (auto& thread : m_workers){
            thread.join();
        }
        m_workers.clear();
    }

    void operator=(const WorkQueue&) = delete;
    WorkQueue(const WorkQueue&) = delete;
};

