function jiggleMouse(pixels, pauseTime)
    if nargin < 1, pixels = 20; end
    if nargin < 2, pauseTime = 0.01; end

    % Get current mouse position (MATLAB: origin bottom-left)
    pos = get(0, 'PointerLocation');

    % Get screen height (Java: origin top-left)
    toolkit = java.awt.Toolkit.getDefaultToolkit();
    screenSize = toolkit.getScreenSize();
    screenHeight = screenSize.getHeight();

    % Convert Y coordinate for Java Robot
    javaY = screenHeight - pos(2);

    robot = java.awt.Robot();

    % Move mouse slightly (e.g., right)
    robot.mouseMove(pos(1) + pixels, javaY);

    pause(pauseTime);

    % Move mouse back
    robot.mouseMove(pos(1), javaY);
end