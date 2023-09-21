//
//  timer.h
//  LE-LP
//
//  Created by moulay abdella chkifa on 9/18/21.
//  Copyright © 2021 moulay abdella chkifa. All rights reserved.
//

#ifndef timer_h
#define timer_h

struct timer {
    
    std::chrono::time_point<std::chrono::steady_clock> start,end;
    std::chrono::duration<float> duration;
    std::string timer_tag;
    
    size_t timer_value(){
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        return 1000*duration.count();
    }
    
    timer(std::string str ):timer_tag{str}{
        start = std::chrono::high_resolution_clock::now();
        std::cout << "Timer " << timer_tag <<  " started: " << std::endl;
    }
    ~timer(){
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Time taken for execution is : " << 1000*duration.count() << " milliseconds" << std::endl;
        std::cout << "Timer " << timer_tag <<  " ended: " << std::endl;
        
    }
};


#endif /* timer_h */
