import time

def run_pubsub_thread(redisdb):
    pubsub = redisdb.pubsub()
    print("run event handler asdasd")
    pubsub.psubscribe(**{'__key*__:*':event_handler}) #**{'__keyevent@0__:expired':self.event_handler}
    thread = pubsub.run_in_thread(sleep_time=1)

def event_handler(msg):
    print("A msg from evt handler!, ",msg)

'''
pubsub = db.pubsub()

print("modelw here")

def event_handler(msg):
    print('Handler', msg)

pubsub.psubscribe(**{'__keyspace@0__:*': event_handler})
'''
