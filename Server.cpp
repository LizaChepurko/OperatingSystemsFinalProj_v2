#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <unordered_map>
#include <atomic>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <functional>
#include <cstring>
#include <sys/wait.h>
#include <fcntl.h>



#define PORT 8080
#define MAX_CLIENTS 30
#define BUFFER_SIZE 2048

/**
 * Active Object for Pipeline Pattern
 * class consists of a task queue, a worker thread, and synchronization primitives to manage task execution.
 * contains its own thread of control and a queue of pending tasks.
 * provides a solid foundation for asynchronously processing tasks in a way
 * that ensures thread safety and efficient resource handling.
 * */
class ActiveObject
{
private:
    std::queue<std::function<void()>> tasks;
    std::mutex mutex;
    std::condition_variable condition;
    std::atomic<bool> running{true}; // Provides atomic operations on particular data types, ensuring thread safety.
    std::thread worker;              // Represents a single thread of execution.

    //This is the function run by the worker thread. 
    //It waits for tasks to be added to the queue and processes them one by one.
    void workerFunction()
    {
        while (running)
        {
            std::function<void()> task;
            {
                std::unique_lock<std::mutex> lock(mutex);
                condition.wait(lock, [this]
                               { return !tasks.empty() || !running; });
                if (!running && tasks.empty())
                    return;
                task = std::move(tasks.front());
                tasks.pop();
            }
            task();
        }
    }

public:
    ActiveObject() : worker(&ActiveObject::workerFunction, this) {}
    //stops the worker thread and waits for it to join, 
    //ensuring that all tasks are completed before the object is destroyed.
    ~ActiveObject()
    {
        {
            std::lock_guard<std::mutex> lock(mutex);
            running = false;
        }
        condition.notify_one();
        worker.join();
    }
    /**
     * Adds a new task to the queue. Notifies the worker thread that a new task is available.
     */
    template <class F>
    void enqueue(F &&f)
    {
        {
            std::lock_guard<std::mutex> lock(mutex);
            tasks.emplace(std::forward<F>(f));
        }
        condition.notify_one();
    }
};

/**
 * Leader-Follower Thread Pool
 * An implementation of the Leader-Follower pattern using a thread pool.
 * The Leader-Follower pattern decouples task submission from task execution using a fixed number of worker threads in a pool.
 * Threads take turns to become the "leader", which picks tasks, while others wait, improving CPU utilization and load balancing.
 */
class LeaderFollowerPool
{
private:
    std::vector<std::thread> threads; // Stores the worker threads.
    std::queue<std::function<void()>> tasks;
    std::mutex mutex;
    std::condition_variable condition; //Notifies that the shared resource is free to access
    bool stop;

public:
    // initializes and starts the worker threads
    LeaderFollowerPool(size_t numThreads) : stop(false)
    {
        for (size_t i = 0; i < numThreads; ++i)
        {
            threads.emplace_back([this]
                                 {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->mutex);
                        this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
                        if (this->stop && this->tasks.empty()) return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    task();
                } });
        }
    }
    //stops the threads and waits for them to join, ensuring clean shutdown
    ~LeaderFollowerPool()
    {
        {
            std::unique_lock<std::mutex> lock(mutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread &worker : threads)
        {
            worker.join();
        }
    }
    //Adds a new task to the queue.
    template <class F>
    void enqueue(F &&f)
    {
        {
            std::unique_lock<std::mutex> lock(mutex);
            tasks.emplace(std::forward<F>(f));
        }
        condition.notify_one();
    }
};

class MSTServer
{
private:
    int graphCommandPipeFd[2]; // Pipe for sending commands to the Graph process
    int graphResultPipeFd[2]; // Pipe for receiving results from the Graph process
    pid_t graphPid;
    ActiveObject pipelineExecutor;
    LeaderFollowerPool lfPool;
    std::mutex graphsMutex;

    void startGraphProcess() {
        // Create command pipe
        if (pipe(graphCommandPipeFd) == -1) {
            perror("pipe (command)");
            exit(EXIT_FAILURE);
        }

        // Create result pipe
        if (pipe(graphResultPipeFd) == -1) {
            perror("pipe (result)");
            exit(EXIT_FAILURE);
        }

        graphPid = fork();
        if (graphPid == -1) {
            perror("fork");
            exit(EXIT_FAILURE);
        } else if (graphPid == 0) { // Child process
            // Close unused ends of the pipes
            close(graphCommandPipeFd[1]);
            close(graphResultPipeFd[0]); 

            // Redirect stdin and stdout to the pipes
            dup2(graphCommandPipeFd[0], STDIN_FILENO);
            dup2(graphResultPipeFd[1], STDOUT_FILENO);
            dup2(graphResultPipeFd[1], STDERR_FILENO);

            // Close the read end of the command pipe
            close(graphCommandPipeFd[0]);
            // Close the write end of the result pipe
            close(graphResultPipeFd[1]);

            execl("./Graph", "Graph", nullptr);
            perror("execl");
            exit(EXIT_FAILURE);
        } else { // Parent process
            // Close unused ends of the pipes
            close(graphCommandPipeFd[0]);
            close(graphResultPipeFd[1]);
        }
    }


        void restartGraphProcess() {
        close(graphResultPipeFd[1]);
        waitpid(graphPid, nullptr, 0);
        startGraphProcess();
    }

public:

    MSTServer(size_t numThreads) : lfPool(numThreads) {
        startGraphProcess();
    }

void handleClient(int clientSocket) {
    char buffer[BUFFER_SIZE] = {0};
    while (true) {
        memset(buffer, 0, BUFFER_SIZE);
        int valread = read(clientSocket, buffer, BUFFER_SIZE);
        if (valread <= 0)
            break;

        std::string request(buffer);
        std::cout << "Received request from client: " << request << std::endl;
        std::istringstream iss(request);
        std::string command;
        iss >> command;

        if (command == "ADD_GRAPH") {
            std::string result = executeCommand(request);
            if (result.empty() || result == "No response from Graph process") {
                result = "Graph added successfully\n";
            }
            std::cout << "Sending to client: " << result << std::endl;
            send(clientSocket, result.c_str(), result.length(), 0);
        }
        else if (command == "SOLVE_MST_PIPELINE" || command == "SOLVE_MST_LF") {
            std::string solveCommand = "SOLVE_MST" + request.substr(command.length());
            
            auto task = [this, clientSocket, solveCommand]() {
                std::string result = executeCommand(solveCommand);
                if (result.empty() || result == "No response from Graph process") {
                    result = "Error: No response from Graph process\n";
                }
                std::cout << "Sending to client: " << result << std::endl;
                send(clientSocket, result.c_str(), result.length(), 0);
            };

            if (command == "SOLVE_MST_PIPELINE") {
                pipelineExecutor.enqueue(task);
            } else {
                lfPool.enqueue(task);
            }
        }
        else {
            std::string errorMsg = "Unknown command\n";
            std::cout << "Sending to client: " << errorMsg << std::endl;
            send(clientSocket, errorMsg.c_str(), errorMsg.length(), 0);
        }
    }
    close(clientSocket);
}


    std::string executeCommand(const std::string& command) {
        std::cout << "Executing command: " << command << std::endl;
        // Write command to the command pipe
        write(graphCommandPipeFd[1], (command + "\n").c_str(), command.length() + 1);

        // Read response from the result pipe
        char buffer[BUFFER_SIZE];
        std::string result;
        ssize_t bytesRead;

        while (true) {
            bytesRead = read(graphResultPipeFd[0], buffer, BUFFER_SIZE - 1);
            if (bytesRead > 0) {
                buffer[bytesRead] = '\0';
                result += buffer;
                // Check if we have reached the end of the response
                if (result.find("\n") != std::string::npos) {
                    break;
                }
            } else if (bytesRead == 0) {
                break; // Reached EOF
            } else {
                perror("read error (result pipe)");
                break;
            }
        }

        std::cout << "Received result from Graph process" << result<< std::endl;
        return result; // No need for success check as the entire response is read
    }


    void run()
    {
        int server_fd, new_socket;
        struct sockaddr_in address;
        int opt = 1;
        int addrlen = sizeof(address);

        if ((server_fd = socket(AF_INET, SOCK_STREAM, 0)) == 0)
        {
            perror("socket failed");
            exit(EXIT_FAILURE);
        }

        if (setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, &opt, sizeof(opt)))
        {
            perror("setsockopt");
            exit(EXIT_FAILURE);
        }

        address.sin_family = AF_INET;
        address.sin_addr.s_addr = INADDR_ANY;
        address.sin_port = htons(PORT);

        if (bind(server_fd, (struct sockaddr *)&address, sizeof(address)) < 0)
        {
            perror("bind failed");
            exit(EXIT_FAILURE);
        }

        if (listen(server_fd, 3) < 0)
        {
            perror("listen");
            exit(EXIT_FAILURE);
        }

        std::cout << "Server is running on port " << PORT << std::endl;

        while (true)
        {
            if ((new_socket = accept(server_fd, (struct sockaddr *)&address, (socklen_t *)&addrlen)) < 0)
            {
                perror("accept");
                continue;
            }

            std::thread clientThread(&MSTServer::handleClient, this, new_socket);
            clientThread.detach();
        }
    }


    ~MSTServer() {
        close(graphCommandPipeFd[1]);
        close(graphResultPipeFd[0]);
        waitpid(graphPid, nullptr, 0);
    }

};

int main()
{
    MSTServer server(5); // Create server with 4 threads in the Leader-Follower pool
    server.run();
    return 0;
}