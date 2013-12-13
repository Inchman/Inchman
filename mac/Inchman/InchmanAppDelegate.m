//
//  InchmanAppDelegate.m
//  Inchman
//
//  Created by Aidan Lane on 23 Jan, 2013
//  Copyright 2011 Monash University. All rights reserved.
//

#import "InchmanAppDelegate.h"

@implementation InchmanAppDelegate


- (void)applicationDidFinishLaunching:(NSNotification *)aNotification {
}

-(void)awakeFromNib {
    [[NSAppleEventManager sharedAppleEventManager] setEventHandler:self andSelector:@selector(handleURLEvent:withReplyEvent:) forEventClass:kInternetEventClass andEventID:kAEGetURL];
    
    // TODO: offer to open an existing xml file and close the dialog immediately if we get a handleURLEvent
}

- (void)handleURLEvent:(NSAppleEventDescriptor*)event withReplyEvent:(NSAppleEventDescriptor*)replyEvent
{
    NSString* url = [[event paramDescriptorForKeyword:keyDirectObject] stringValue];
    NSLog(@"%@", url);

    NSString *pythonHome = [NSString pathWithComponents:[NSArray arrayWithObjects:[[NSBundle mainBundle] bundlePath], @"Contents", @"Frameworks", @"Python.framework", @"Versions", @"2.7", nil]];
    NSString *pythonExec = [NSString pathWithComponents:[NSArray arrayWithObjects:[[NSBundle mainBundle] bundlePath], @"Contents", @"Frameworks", @"Python.framework", @"Versions", @"2.7", @"bin", @"python", nil]];
    NSString *protocolhandlerScript = [NSString pathWithComponents:[NSArray arrayWithObjects:[[NSBundle mainBundle] bundlePath], @"Contents", @"MacOS", @"protocolhandler.py", nil]];
    NSString *inchmanExec = [NSString pathWithComponents:[NSArray arrayWithObjects:[[NSBundle mainBundle] bundlePath], @"Contents", @"MacOS", @"InchmanExec", nil]];

    NSString *scriptStr = [NSString stringWithFormat:
                           @"set workingPath to (choose folder with prompt \"Choose working directory\")\n"
                           "set posixWorkingPath to POSIX path of workingPath\n"
                           "tell application \"Terminal\"\n"
                           "  activate\n"
                           // FUTURE: ensure that we're using '/bin/sh', as we use it's 'export' command
                           "  do script \"cd '\" & posixWorkingPath & \"'; export PYTHONHOME='%@'; %@ %@ '%@' '%@'\"\n"
                           "end tell", pythonHome, pythonExec, protocolhandlerScript, inchmanExec, url];
    
    NSAppleScript *script = [[NSAppleScript alloc] initWithSource:scriptStr];
    [script executeAndReturnError:nil];
    
    [NSApp terminate:self];
}


@end
