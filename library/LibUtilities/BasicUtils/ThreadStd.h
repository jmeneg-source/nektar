///////////////////////////////////////////////////////////////////////////////
//
// File: ThreadStd.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Thread manager implementation using std::thread.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBUTILITIES_THREADSTD_H_
#define NEKTAR_LIBUTILITIES_THREADSTD_H_

#include <condition_variable>
#include <map>
#include <mutex>
#include <queue>
#include <vector>

#include <LibUtilities/BasicUtils/Thread.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar::Thread
{

typedef std::unique_lock<std::mutex> Lock;

class ThreadWorkerStd;

/**
 * @brief Implementation of ThreadManager using std::thread.
 */
class ThreadManagerStd : public ThreadManager
{
    /// Lightweight barrier class
    class Barrier
    {
    public:
        explicit Barrier(std::size_t iCount)
            : mThreshold(iCount), mCount(iCount), mGeneration(0)
        {
        }

        void Wait()
        {
            std::unique_lock<std::mutex> lLock{mMutex};
            auto lGen = mGeneration;
            if (!--mCount)
            {
                mGeneration++;
                mCount = mThreshold;
                mCond.notify_all();
            }
            else
            {
                mCond.wait(lLock, [this, lGen] { return lGen != mGeneration; });
            }
        }

    private:
        std::mutex mMutex;
        std::condition_variable mCond;
        std::size_t mThreshold;
        std::size_t mCount;
        std::size_t mGeneration;
    };

    /**
     * So the workers can access the master queue and locks.
     */
    friend class ThreadWorkerStd;

public:
    /// Constructs a ThreadManagerStd.
    ThreadManagerStd(unsigned int numWorkers);
    /// Shuts down threading.
    ~ThreadManagerStd() override;

    /// Called by the factory method.
    static ThreadManagerSharedPtr Create(unsigned int numT)
    {
        return std::shared_ptr<ThreadManager>(new ThreadManagerStd(numT));
    }

protected:
    void v_QueueJobs(std::vector<ThreadJob *> &joblist) override;
    void v_QueueJob(ThreadJob *job) override;
    unsigned int v_GetNumWorkers() override;
    unsigned int v_GetWorkerNum() override;
    void v_SetNumWorkers(const unsigned int num) override;
    void v_SetNumWorkers() override;
    unsigned int v_GetMaxNumWorkers() override;
    void v_Wait() override;
    void v_SetChunkSize(unsigned int chnk) override;
    void v_SetSchedType(SchedType s) override;
    bool v_InThread() override;
    void v_Hold() override;
    const std::string &v_GetType() const override;

private:
    ThreadManagerStd();
    ThreadManagerStd(const ThreadManagerStd &);
    bool IsWorking();
    void SetNumWorkersImpl(const unsigned int num);

    // Member variables
    const unsigned int m_numThreads;
    unsigned int m_numWorkers;
    std::queue<ThreadJob *> m_masterQueue;
    std::mutex m_masterQueueMutex;
    std::mutex m_masterActiveMutex;
    std::condition_variable m_masterQueueCondVar;
    std::condition_variable m_masterActiveCondVar;
    ThreadWorkerStd **m_threadList;
    std::thread **m_threadThreadList;
    std::thread::id m_masterThreadId;
    bool *m_threadBusyList;
    bool *m_threadActiveList;
    unsigned int m_chunkSize;
    SchedType m_schedType;
    Barrier *m_barrier;
    std::map<std::thread::id, unsigned int> m_threadMap;
    static std::string className;
    std::string m_type;
};

/**
 * @brief Implementation class for ThreadManagerStd.
 *
 * Each instance of this class corresponds to a worker thread.
 * Instances manage their own queue of jobs to run, grabbing new
 * jobs from the master queue when it is exhausted.
 */
class ThreadWorkerStd
{
public:
    /// Constructor
    ThreadWorkerStd(ThreadManagerStd *threadManager, unsigned int workerNum);
    /// Destructor.
    ~ThreadWorkerStd();
    /// This provides the interface that std::thread uses to start the worker.
    void operator()()
    {
        MainLoop();
    };
    /**
     * @brief Return the index of the worker thread.
     * @return Index of worker thread, an integer between 0 and
     *        (number_of_threads - 1)
     */
    unsigned int GetWorkerNum()
    {
        return m_threadNum;
    };
    /**
     * @brief A signal to shut down.
     *
     * If this method is called the worker will shut down.  Used by the
     * ThreadManagerStd to stop threading.
     */
    void Stop()
    {
        m_keepgoing = false;
    };

private:
    ThreadWorkerStd();
    ThreadWorkerStd(const ThreadWorkerStd &);
    void MainLoop();
    void LoadJobs();
    unsigned int GetNumToLoad();
    void WaitForActive();
    void RunJobs();

    // Member variables
    ThreadManagerStd *m_threadManager;
    std::queue<ThreadJob *> m_workerQueue;
    bool m_keepgoing;
    unsigned int m_threadNum;
};

} // namespace Nektar::Thread

#endif
