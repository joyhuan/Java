/** Tutorial from 1BestCsharp blog
* (https://1bestcsharp.blogspot.com/2017/08/java-login-and-register-form-design.html)
* This is a much nice and more intuitive way to create GUI. When I was 
* a freshman, we still needed to manually set the location of the components.
* Unnecessary toil. */

public class JFrameForm extends javax.swing.JFram{
    public RegisterForm(){
        initComponents(); 
    }

    // jlabelClose = close form on jlabel click
    private void jLabelCloseMouseClicked(java.awt.event.MouseEvent evt) {                                         
        System.exit(0);
    }   

    // jlabelMin = minimize form on jlabel click
    private void jLabelMinMouseClicked(java.awt.event.MouseEvent evt) {                                       
        this.setState(JFrame.ICONIFIED);
    }

    // jLabelRegister = open login form on jlabel click (on register form)
    private void jLabelRegisterMouseClicked(java.awt.event.MouseEvent evt) {                                            
        LoginForm lgf = new LoginForm();
        lgf.setVisible(true);
        lgf.pack();
        lgf.setLocationRelativeTo(null);
        lgf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        this.dispose();
    }

    // // jLabelRegister = open login form on jlabel click  (on login form)
    // private void jLabelRegisterMouseClicked(java.awt.event.MouseEvent evt) {                                            
    //     RegisterForm rgf = new RegisterForm();
    //     rgf.setVisible(true);
    //     rgf.pack();
    //     rgf.setLocationRelativeTo(null);
    //     rgf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    // this.dispose();
    // }    
}
