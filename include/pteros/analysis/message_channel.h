
/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 *
 * Please note, that Artistic License is slightly more restrictive
 * then GPL license in terms of distributing the modified versions
 * of this software (they should be approved first).
 * Read http://www.opensource.org/licenses/artistic-license-2.0.php
 * for details. Such license fits scientific software better then
 * GPL because it prevents the distribution of bugged derivatives.
 *
*/

#ifndef MESSAGE_CHANNEL_H
#define MESSAGE_CHANNEL_H

#include <thread>
#include <queue>

template<class T>
class Message_channel {
public:
    Message_channel(){
        buffer_size = 10;
        stop_requested = false;
    }

    Message_channel(int sz){
        buffer_size = sz;
        stop_requested = false;
    }

    void send_stop(){
        std::lock_guard<std::mutex> lock(mutex);
        stop_requested = true;
        cond.notify_all();
    }

    bool empty(){
        std::lock_guard<std::mutex> lock(mutex);
        return queue.empty();
    }

    void send(T const& data){
        std::unique_lock<std::mutex> lock(mutex);
        // Wait until buffer will clear a bit or until stop is requested
        while(queue.size()>=buffer_size && !stop_requested){
            cond.wait(lock);
        }
        // If stop requested just do nothing
        if(stop_requested) return;

        queue.push(data);
        cond.notify_all();
    }

    bool recieve(T& popped_value){
        std::unique_lock<std::mutex> lock(mutex);
        // Wait until something appears in the queue or until stop requested
        while(queue.empty() && !stop_requested){
            cond.wait(lock);
        }
        // If stop requested see if there is something in, if not return false
        if(stop_requested && queue.empty()) return false;

        popped_value=queue.front();
        queue.pop();
        cond.notify_all();
        return true;
    }

private:
    int buffer_size;
    std::condition_variable cond;
    std::mutex mutex;
    std::queue<T> queue;
    bool stop_requested;
};


#endif // MESSAGE_CHANNEL_H
