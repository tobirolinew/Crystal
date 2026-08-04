      select type(x)
      type is(integer(i32)); rflat=reshape(dble(x),[m,n])
      type is(integer(i64)); rflat=reshape(dble(x),[m,n])
      type is(real(r32)); rflat=reshape(dble(x),[m,n])
      type is(real(r64)); rflat=reshape(dble(x),[m,n])
      class default
         call catch_error(&
                  err=.true.,&
                  msg='unknown array type.',&
                  proc='get_num_lex_sort_perm')
      end select
